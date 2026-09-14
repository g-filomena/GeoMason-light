package sim.graph;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.Stack;
import org.javatuples.Pair;
import org.locationtech.jts.geom.Coordinate;
import org.locationtech.jts.geom.Envelope;
import org.locationtech.jts.index.strtree.ItemDistance;
import org.locationtech.jts.index.strtree.STRtree;
import sim.routing.Astar;
import sim.routing.Route;

public class Islands {

  Set<NodeGraph> visitedNodes = new HashSet<>();
  private Graph graph;
  List<Set<NodeGraph>> islands = new ArrayList<>();

  /**
   * Constructs an Islands object for finding and merging disconnected islands in the given graph.
   *
   * @param graph the graph to analyze for islands
   */
  public Islands(Graph graph) {
    this.graph = graph;
  }

  /**
   * Finds and returns a list of sets, where each set contains nodes that form a disconnected island
   * in the graph.
   *
   * @param edges the set of edges to consider when finding disconnected islands
   * @return a list of sets, each representing a disconnected island of nodes
   */
  public List<Set<NodeGraph>> findDisconnectedIslands(Set<EdgeGraph> edges) {
    // LinkedHashSet: NodeGraph has no hashCode of its own, so a HashSet iterates in identity-hash
    // order, which differs between JVM builds. The search below is seeded from these nodes in turn,
    // so their order becomes the island order, and mergeConnectedIslands bridges the first
    // cross-island edge it finds in it.
    Set<NodeGraph> nodes = new LinkedHashSet<>(GraphUtils.nodesFromEdges(edges));
    this.visitedNodes = new HashSet<>();
    islands.clear();

    // Sequential, deliberately. This was a parallelStream whose entire body sat inside
    // synchronized (visitedNodes), so every DFS serialised anyway and the only thing the
    // fork-join pool added was overhead - on a hot path called once per agent, and again on
    // every iteration of the merge loop below.
    for (NodeGraph currentNode : nodes) {
      if (visitedNodes.contains(currentNode)) {
        continue;
      }
      // Insertion-ordered too: findConnectingBridge stops at the first node adjacent to another
      // island, and the closest-pair search breaks ties on whichever it reaches first.
      Set<NodeGraph> currentIsland = new LinkedHashSet<>();
      dfs(currentNode, currentIsland, nodes, edges);
      islands.add(currentIsland);
    }

    return islands;
  }

  /**
   * Performs a depth-first search (DFS) to explore and mark all nodes connected to the given node.
   *
   * @param node the starting node for the DFS
   * @param currentIsland the set to store nodes that are part of the current island
   * @param nodes the set of all nodes to consider in the DFS
   * @param edges the set of edges to consider in the DFS
   */
  private void dfs(NodeGraph node, Set<NodeGraph> currentIsland, Set<NodeGraph> nodes,
      Set<EdgeGraph> edges) {

    Stack<NodeGraph> stack = new Stack<>();
    stack.push(node);

    while (!stack.isEmpty()) {
      NodeGraph currentNode = stack.pop();

      if (!visitedNodes.contains(currentNode)) {
        visitedNodes.add(currentNode);
        currentIsland.add(currentNode);
        for (NodeGraph adjacentNode : currentNode.getAdjacentNodes()) {
          if (!visitedNodes.contains(adjacentNode)
              && edges.contains(graph.getEdgeBetween(currentNode, adjacentNode))) {
            stack.push(adjacentNode);
          }
        }
      }
    }
  }

  /**
   * Merges disconnected islands in the graph by adding edges to connect them, ensuring all nodes
   * are connected.
   *
   * @param edges the set of edges to consider when merging islands
   * @return the updated set of edges with added bridges to connect islands
   */
  public Set<EdgeGraph> mergeConnectedIslands(Set<EdgeGraph> edges) {

    // Insertion-ordered: these go into the caller's set and are iterated again downstream.
    Set<EdgeGraph> bridges = new LinkedHashSet<>();
    islands = findDisconnectedIslands(edges);
    if (islands.size() == 1) {
      return edges;
    }
    Astar aStar = new Astar();
    // Bounded. If A* cannot route between the closest pair - the two islands sit in parts of the
    // graph with no path between them - the bridge adds no edges, the island count does not fall,
    // and this loop never ends. It is a hang rather than a crash, which is worse: the run simply
    // stops making progress. Each pass must strictly reduce the island count or we stop and hand
    // back what we have, which is a usable if imperfectly connected known network.
    int islandsBefore = islands.size();
    int stalls = 0;
    while (islands.size() > 1 && stalls < 3) {
      EdgeGraph connectingEdge = findConnectingBridge(islands, graph);
      if (connectingEdge != null) {
        bridges.add(connectingEdge);
      } else {
        Pair<NodeGraph, NodeGraph> closestNodes = findClosestPairAcrossAllIslands();
        if (closestNodes == null) {
          break;
        }
        Route route =
            aStar.astarRoute(closestNodes.getValue0(), closestNodes.getValue1(), graph, null);
        if (route == null || route.edgesSequence == null || route.edgesSequence.isEmpty()) {
          stalls++;
          continue;
        }
        bridges.addAll(route.edgesSequence);
      }
      edges.addAll(bridges);
      islands = findDisconnectedIslands(edges);
      if (islands.size() >= islandsBefore) {
        stalls++;
      } else {
        stalls = 0;
        islandsBefore = islands.size();
      }
    }
    if (islands.size() > 1) {
      // Counted, not silent: a known network left in pieces is a network some origin-destination
      // pairs have no route through, and that surfaces far away as a route-choice model quietly
      // substituting another one.
      incompleteMerges++;
    }
    return edges;
  }

  /** Merges that gave up with the islands still separate. */
  private static volatile long incompleteMerges = 0;

  /** How many merges have given up so far, across the run. */
  public static long incompleteMerges() {
    return incompleteMerges;
  }

  /**
   * Finds the globally closest pair of nodes lying in two different islands.
   *
   * <p>Plain nested loops over each unordered island pair. What this replaces was a nested
   * {@code parallelStream} that materialised a {@code Pair} and a {@code SimpleEntry} for every
   * cross-island node pair before taking the minimum, and enumerated each unordered pair twice
   * because it never excluded the symmetric case. On the activity module - whose agents anchor on
   * home, work and two persona destinations, so their known space fragments into several
   * well-separated islands - that allocation was the single hottest thing in the model, and it is
   * paid once per agent.
   *
   *
   * @return the closest cross-island pair, or null when there are fewer than two islands
   */
  private Pair<NodeGraph, NodeGraph> findClosestPairAcrossAllIslands() {
    if (islands.size() < 2) {
      return null;
    }

    // One spatial index per island, then for each unordered pair walk the smaller island and ask
    // the larger one's index for its nearest node. Closest-pair-between-two-point-sets is what a
    // spatial index is for; the nested loops this replaces compared every node of every island
    // against every node of every other, which is quadratic in the size of an agent's known space
    // and is paid once per agent.
    List<STRtree> trees = new ArrayList<>(islands.size());
    for (Set<NodeGraph> island : islands) {
      STRtree tree = new STRtree();
      for (NodeGraph node : island) {
        Coordinate coordinate = node.getCoordinate();
        tree.insert(new Envelope(coordinate), node);
      }
      tree.build();
      trees.add(tree);
    }

    ItemDistance distance =
        (a, b) -> GraphUtils.nodesDistance((NodeGraph) a.getItem(), (NodeGraph) b.getItem());

    NodeGraph bestFrom = null;
    NodeGraph bestTo = null;
    double bestDistance = Double.MAX_VALUE;

    for (int i = 0; i < islands.size(); i++) {
      for (int j = i + 1; j < islands.size(); j++) {
        int from = islands.get(i).size() <= islands.get(j).size() ? i : j;
        int to = from == i ? j : i;
        STRtree target = trees.get(to);
        for (NodeGraph node : islands.get(from)) {
          Object nearest =
              target.nearestNeighbour(new Envelope(node.getCoordinate()), node, distance);
          if (nearest == null) {
            continue;
          }
          NodeGraph other = (NodeGraph) nearest;
          double d = GraphUtils.nodesDistance(node, other);
          if (d < bestDistance) {
            bestDistance = d;
            bestFrom = node;
            bestTo = other;
          }
        }
      }
    }
    return bestFrom == null ? null : new Pair<>(bestFrom, bestTo);
  }

  /**
   * Finds an edge that can act as a bridge to connect two disconnected islands in the graph.
   *
   * <p>Considers each unordered island pair once, and stops at the first edge it finds.
   *
   * @param islands the list of sets of nodes representing the disconnected islands
   * @param graph the graph to analyze for potential connecting bridges
   * @return an edge that connects two islands, or null if no such edge exists
   */
  private static EdgeGraph findConnectingBridge(List<Set<NodeGraph>> islands, Graph graph) {
    // An edge between two nodes exists only if they are adjacent, so only adjacent pairs are worth
    // asking about. Walking each node's own neighbours is O(N x degree), where a street node has
    // three or four, rather than the O(N^2) of every node in every island against every node in
    // every other.
    Map<NodeGraph, Integer> islandOf = new HashMap<>();
    for (int i = 0; i < islands.size(); i++) {
      for (NodeGraph node : islands.get(i)) {
        islandOf.put(node, i);
      }
    }
    for (int i = 0; i < islands.size(); i++) {
      for (NodeGraph node : islands.get(i)) {
        for (NodeGraph adjacent : node.getAdjacentNodes()) {
          Integer other = islandOf.get(adjacent);
          if (other != null && other != i) {
            EdgeGraph edge = graph.getEdgeBetween(node, adjacent);
            if (edge != null) {
              return edge;
            }
          }
        }
      }
    }
    return null;
  }
}
