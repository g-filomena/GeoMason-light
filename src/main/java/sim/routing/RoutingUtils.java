package sim.routing;

import java.util.ArrayList;
import java.util.List;
import java.util.stream.Collectors;
import org.locationtech.jts.planargraph.DirectedEdge;
import sim.graph.EdgeGraph;
import sim.graph.NodeGraph;

public class RoutingUtils {

  /**
   * Identifies the previous junction traversed in a dual graph path to avoid traversing an
   * unnecessary segment in the primal graph.
   *
   * @param sequenceDirectedEdges A sequence of GeomPlanarGraphDirectedEdge representing the path.
   * @return The previous junction node.
   */
  public static NodeGraph getPreviousJunction(List<DirectedEdge> sequenceDirectedEdges) {

    if (sequenceDirectedEdges.size() == 1) {
      return (NodeGraph) sequenceDirectedEdges.get(0).getFromNode();
    }

    int ixLast = sequenceDirectedEdges.size() - 1;
    int ixBeforeLast = sequenceDirectedEdges.size() - 2;
    NodeGraph lastCentroid =
        ((EdgeGraph) sequenceDirectedEdges.get(ixLast).getEdge()).getDualNode();
    NodeGraph otherCentroid =
        ((EdgeGraph) sequenceDirectedEdges.get(ixBeforeLast).getEdge()).getDualNode();
    return getPrimalJunction(lastCentroid, otherCentroid);
  }

  /**
   * Given two centroids (nodes in the dual graph), identifies their shared junction (i.e., the
   * junction shared by the corresponding primal segments).
   *
   * @param centroid A dual node.
   * @param otherCentroid Another dual node.
   * @return The common primal junction node.
   */
  public static NodeGraph getPrimalJunction(NodeGraph centroid, NodeGraph otherCentroid) {

    EdgeGraph edge = centroid.getPrimalEdge();
    EdgeGraph otherEdge = otherCentroid.getPrimalEdge();

    if (edge.getFromNode().equals(otherEdge.getFromNode())
        || edge.getFromNode().equals(otherEdge.getToNode())) {
      return edge.getFromNode();
    } else if (edge.getToNode().equals(otherEdge.getFromNode())
        || edge.getToNode().equals(otherEdge.getToNode())) {
      return edge.getToNode();
    } else {
      return null;
    }
  }

  /**
   * Returns all the primal nodes traversed in a path.
   *
   * @param directedEdgesSequence A sequence of GeomPlanarGraphDirectedEdge representing the path.
   * @return A list of primal nodes.
   */
  public static List<NodeGraph> getNodesFromDirectedEdgesSequence(
      List<DirectedEdge> directedEdgesSequence) {

    List<NodeGraph> nodesSequence = new ArrayList<>();
    if (directedEdgesSequence.isEmpty()) {
      return nodesSequence;
    }

    nodesSequence = directedEdgesSequence.stream()
        .map(directedEdge -> (NodeGraph) directedEdge.getFromNode()).collect(Collectors.toList());

    DirectedEdge lastEdge = directedEdgesSequence.get(directedEdgesSequence.size() - 1);
    nodesSequence.add((NodeGraph) lastEdge.getToNode());
    return nodesSequence;
  }

  /**
   * Returns all the centroids (nodes in the dual graph) traversed in a path.
   *
   * @param sequenceDirectedEdges A sequence of GeomPlanarGraphDirectedEdge representing the path.
   * @return A list of centroids (dual nodes).
   */
  public static List<NodeGraph> getCentroidsFromEdgesSequence(
      List<DirectedEdge> sequenceDirectedEdges) {
    return sequenceDirectedEdges.stream()
        .map(planarDirectedEdge -> ((EdgeGraph) planarDirectedEdge.getEdge()).getDualNode())
        .collect(Collectors.toList());
  }

}
