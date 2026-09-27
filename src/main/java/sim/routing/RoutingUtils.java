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

    // Walk the sequence junction by junction: where two consecutive segments are parallel (they
    // share both ends), only the junction the walk arrived by tells which end it crossed at.
    NodeGraph junction = null;
    for (int i = 1; i < sequenceDirectedEdges.size(); i++) {
      NodeGraph centroid = ((EdgeGraph) sequenceDirectedEdges.get(i - 1).getEdge()).getDualNode();
      NodeGraph nextCentroid = ((EdgeGraph) sequenceDirectedEdges.get(i).getEdge()).getDualNode();
      junction = getPrimalJunction(centroid, nextCentroid, junction);
    }
    return junction;
  }

  /**
   * Identifies the junction at which a walk moves from one segment to the next, given the junction
   * by which it arrived on the first segment. The walk leaves a segment at its far end, so that end
   * is answered when the next segment shares it. This is what tells parallel segments (which share
   * both ends) apart; for any other pair it gives what
   * {@link #getPrimalJunction(NodeGraph, NodeGraph)} gives.
   *
   * @param centroid The dual node of the segment the walk is on.
   * @param otherCentroid The dual node of the next segment.
   * @param arrivalJunction The junction by which the walk arrived on {@code centroid}'s segment;
   *        null when unknown.
   * @return The junction between the two segments, or null if they share none.
   */
  public static NodeGraph getPrimalJunction(NodeGraph centroid, NodeGraph otherCentroid,
      NodeGraph arrivalJunction) {

    EdgeGraph edge = centroid.getPrimalEdge();
    EdgeGraph otherEdge = otherCentroid.getPrimalEdge();
    NodeGraph farEnd = null;
    if (edge.getFromNode().equals(arrivalJunction)) {
      farEnd = edge.getToNode();
    } else if (edge.getToNode().equals(arrivalJunction)) {
      farEnd = edge.getFromNode();
    }
    if (farEnd != null
        && (farEnd.equals(otherEdge.getFromNode()) || farEnd.equals(otherEdge.getToNode()))) {
      return farEnd;
    }
    return getPrimalJunction(centroid, otherCentroid);
  }

  /**
   * Given two centroids (nodes in the dual graph), identifies their shared junction (i.e., the
   * junction shared by the corresponding primal segments). Parallel segments share both ends, and
   * this answers the first segment's from-node; where the walk's direction is known, use
   * {@link #getPrimalJunction(NodeGraph, NodeGraph, NodeGraph)}.
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
