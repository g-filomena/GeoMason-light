package sim.routing;

import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.PriorityQueue;
import java.util.Set;
import java.util.function.Predicate;
import org.locationtech.jts.geom.Coordinate;
import org.locationtech.jts.planargraph.DirectedEdge;
import sim.graph.EdgeGraph;
import sim.graph.Graph;
import sim.graph.NodeGraph;
import sim.util.geo.GeometryUtilities;

/**
 * A* over a {@link Graph}, minimising metric length with a straight-line heuristic.
 *
 * <p>Edges may be excluded, either as a set of edge IDs or as a {@link Predicate} tested on the
 * edge. The predicate form is the cheaper one when the exclusion is derived rather than
 * enumerated: only the edges the search actually expands are ever tested, whereas a set has to be
 * built for the whole graph first.
 *
 * <p>Searches may have several targets. {@link #astarRouteAllowing(NodeGraph, Collection, Graph,
 * Predicate)} returns the route to whichever target is settled first, which is one search rather
 * than one per candidate.
 *
 * <p>The predicate methods carry their own name rather than overloading {@code astarRoute}: a
 * caller passing a literal {@code null} could not tell the two apart.
 */
public class Astar {

  /** The target settled by the last multi-target search, or null. */
  private NodeGraph reachedTarget;

  /**
   * Finds the shortest route between two nodes, avoiding a set of edge IDs.
   *
   * @param originNode the starting node for the route
   * @param destinationNode the destination node for the route
   * @param graph the graph on which to calculate the route
   * @param edgesToAvoid IDs of edges the route may not use; null or empty avoids nothing
   * @return the route, or null if no path exists
   */
  public Route astarRoute(NodeGraph originNode, NodeGraph destinationNode, Graph graph,
      Set<Integer> edgesToAvoid) {

    final Predicate<EdgeGraph> allowed = (edgesToAvoid == null || edgesToAvoid.isEmpty())
        ? edge -> true
        : edge -> !edgesToAvoid.contains(edge.getID());
    return astarRouteAllowing(originNode, destinationNode, graph, allowed);
  }

  /**
   * Finds the shortest route between two nodes, admitting only edges the predicate accepts.
   *
   * @param originNode the starting node for the route
   * @param destinationNode the destination node for the route
   * @param graph the graph on which to calculate the route
   * @param edgeAllowed tested on each edge as the search reaches it
   * @return the route, or null if no path exists
   */
  public Route astarRouteAllowing(NodeGraph originNode, NodeGraph destinationNode, Graph graph,
      Predicate<EdgeGraph> edgeAllowed) {
    return astarRouteAllowing(originNode, Collections.singleton(destinationNode), graph,
        edgeAllowed);
  }

  /**
   * Finds the shortest route from {@code originNode} to whichever of {@code targets} is settled
   * first.
   *
   * <p>One search, not one per target: A* settles nodes in order of {@code g + h}, so the first
   * target polled is the one the search reaches first under the heuristic. {@link #reachedTarget()}
   * names it afterwards. With several targets the heuristic is the distance to the nearest of them,
   * which keeps it admissible.
   *
   * @param originNode the starting node for the route
   * @param targets the acceptable end nodes; the search stops at the first one settled
   * @param graph the graph on which to calculate the route
   * @param edgeAllowed tested on each edge as the search reaches it
   * @return the route to the first target reached, or null if none is reachable
   */
  public Route astarRouteAllowing(NodeGraph originNode, Collection<NodeGraph> targets, Graph graph,
      Predicate<EdgeGraph> edgeAllowed) {

    reachedTarget = null;
    if (originNode == null || targets == null || targets.isEmpty()) {
      return null;
    }
    final Set<NodeGraph> targetSet =
        (targets instanceof Set) ? (Set<NodeGraph>) targets : new HashSet<>(targets);

    final Map<NodeGraph, NodeWrapper> nodeWrappersMap = new HashMap<>();
    final PriorityQueue<Entry> openSet =
        new PriorityQueue<>(Comparator.comparingDouble(entry -> entry.fx));
    final Set<NodeGraph> closedSet = new HashSet<>();
    // NodeWrapper.hx defaults to 0, which is indistinguishable from a computed zero, so the
    // heuristic is remembered here rather than inferred from the wrapper.
    final Map<NodeGraph, Double> heuristics = new HashMap<>();

    final NodeWrapper originWrapper = getNodeWrapper(originNode, nodeWrappersMap);
    originWrapper.gx = 0.0;
    originWrapper.hx = heuristic(originNode, targetSet, heuristics);
    originWrapper.fx = originWrapper.hx;
    openSet.add(new Entry(originNode, originWrapper.fx));

    while (!openSet.isEmpty()) {

      final Entry entry = openSet.poll();
      final NodeGraph currentNode = entry.node;

      // Lazy deletion: a node may sit in the queue several times, once per improvement. The first
      // time it is polled it is final, and later entries are stale. Testing membership in the queue
      // instead - and removing from it - is a linear scan of the queue per improved neighbour.
      if (!closedSet.add(currentNode)) {
        continue;
      }

      if (targetSet.contains(currentNode)) {
        reachedTarget = currentNode;
        return reconstructPath(nodeWrappersMap.get(currentNode));
      }

      final NodeWrapper currentWrapper = nodeWrappersMap.get(currentNode);

      // The node's own outgoing directed edges carry the target node, the undirected edge and the
      // directed edge together, so no getEdgeBetween/getDirectedEdgeBetween lookups are needed -
      // each of those allocates a key pair, and the old shape performed three per neighbour.
      for (DirectedEdge outEdge : currentNode.getOutDirectedEdges()) {
        final NodeGraph targetNode = (NodeGraph) outEdge.getToNode();
        if (closedSet.contains(targetNode)) {
          continue;
        }
        final EdgeGraph edge = (EdgeGraph) outEdge.getEdge();
        if (!edgeAllowed.test(edge)) {
          continue;
        }

        final double tentativeGx = currentWrapper.gx + edge.getLength();
        final NodeWrapper nextWrapper = getNodeWrapper(targetNode, nodeWrappersMap);
        if (tentativeGx >= nextWrapper.gx) {
          continue;
        }

        nextWrapper.previousWrapper = currentWrapper;
        nextWrapper.directedEdgeFrom = outEdge;
        nextWrapper.gx = tentativeGx;
        nextWrapper.hx = heuristic(targetNode, targetSet, heuristics);
        nextWrapper.fx = tentativeGx + nextWrapper.hx;
        openSet.add(new Entry(targetNode, nextWrapper.fx));
      }
    }
    return null;
  }

  /**
   * The target settled by the last multi-target search.
   *
   * @return the node reached, or null if the last search found no path
   */
  public NodeGraph reachedTarget() {
    return reachedTarget;
  }

  /** A queue entry pairing a node with the cost it was enqueued at, so entries can go stale. */
  private static final class Entry {
    private final NodeGraph node;
    private final double fx;

    private Entry(NodeGraph node, double fx) {
      this.node = node;
      this.fx = fx;
    }
  }

  private static NodeWrapper getNodeWrapper(NodeGraph node,
      Map<NodeGraph, NodeWrapper> nodeWrappersMap) {
    return nodeWrappersMap.computeIfAbsent(node, NodeWrapper::new);
  }

  /**
   * Rebuilds the path by walking the predecessors back from the end node.
   *
   * @param endWrapper the wrapper at the node the search settled on
   * @return the route from the origin to that node
   */
  private static Route reconstructPath(NodeWrapper endWrapper) {

    final Route route = new Route();
    final List<DirectedEdge> directedEdgesSequence = new ArrayList<>();

    // Appended then reversed once: inserting at the head of an ArrayList shifts the whole list on
    // every step, which is quadratic in the length of the route.
    NodeWrapper wrapper = endWrapper;
    while (wrapper != null && wrapper.previousWrapper != null) {
      directedEdgesSequence.add(wrapper.directedEdgeFrom);
      wrapper = wrapper.previousWrapper;
    }
    Collections.reverse(directedEdgesSequence);

    route.directedEdgesSequence = directedEdgesSequence;
    if (!route.directedEdgesSequence.isEmpty()) {
      route.computeRouteSequences();
    }
    return route;
  }

  /**
   * Straight-line distance to the nearest of the targets, which keeps the heuristic admissible when
   * there is more than one.
   *
   * @param node the node to score
   * @param targets the acceptable end nodes
   * @param cache distances already computed in this search
   * @return the distance to the closest target
   */
  private static double heuristic(NodeGraph node, Set<NodeGraph> targets,
      Map<NodeGraph, Double> cache) {

    final Double remembered = cache.get(node);
    if (remembered != null) {
      return remembered;
    }
    final Coordinate nodeCoords = node.getCoordinate();
    double best = Double.MAX_VALUE;
    for (NodeGraph target : targets) {
      final double distance =
          GeometryUtilities.euclideanDistance(nodeCoords, target.getCoordinate());
      if (distance < best) {
        best = distance;
      }
    }
    cache.put(node, best);
    return best;
  }
}
