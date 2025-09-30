package sim.routing;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;
import org.locationtech.jts.geom.Coordinate;
import org.locationtech.jts.geom.GeometryFactory;
import org.locationtech.jts.geom.LineString;
import org.locationtech.jts.planargraph.DirectedEdge;
import sim.graph.EdgeGraph;
import sim.graph.Graph;
import sim.graph.NodeGraph;

/**
 * A class for storing the sequence of GeomPlanarGraphDirectedEdge in a path and the sequence of
 * NodeWrappers. It supports path algorithms and provides some utilities.
 */
public class Route {

  // always primal
  public NodeGraph originNode;
  public NodeGraph destinationNode;
  public List<DirectedEdge> directedEdgesSequence = new ArrayList<>();
  public List<NodeGraph> nodesSequence = new ArrayList<>();
  public List<EdgeGraph> edgesSequence = new ArrayList<>();
  public List<NodeGraph> dualNodesSequence = new ArrayList<>();
  private LineString lineString;

  GeometryFactory FACTORY = new GeometryFactory();
  public Map<String, Object> attributes = new HashMap<>();
  private Set<NodeGraph> visitedLocations = new HashSet<>();
  public boolean social = false;
  private double length;

  public Route() {

  }

  /**
   * Constructs a Route object with a given sequence of directed edges.
   *
   * @param directedEdgeSequence A list of DirectedEdge objects representing the path.
   */
  public Route(List<DirectedEdge> directedEdgeSequence) {
    this.directedEdgesSequence = directedEdgeSequence;
  }

  /**
   * Resets the route with a new sequence of directed edges and recomputes the route sequences.
   *
   * @param directedEdgeSequence A list of DirectedEdge objects representing the new path.
   */
  public void resetRoute(List<DirectedEdge> directedEdgeSequence) {
    directedEdgesSequence = new ArrayList<>(directedEdgeSequence);
    nodesSequence = new ArrayList<>();
    edgesSequence = new ArrayList<>();
    computeRouteSequences();
  }

  /**
   * Computes the sequences of nodes and edges for the route, setting the origin and destination
   * nodes, and obtaining the route geometry.
   */
  public void computeRouteSequences() {

    nodesSequence();
    edgeSequence();
    originNode = nodesSequence.get(0);
    destinationNode = nodesSequence.get(nodesSequence.size() - 1);
    obtainRouteLineGeometry();
  }

  /**
   * Returns all the primal nodes traversed in a path.
   *
   * @param directedEdgesSequence A sequence of GeomPlanarGraphDirectedEdge representing the path.
   * @return A list of primal nodes.
   */
  private void nodesSequence() {
    if (directedEdgesSequence.isEmpty()) {
      return;
    }

    nodesSequence = directedEdgesSequence.stream()
        .map(directedEdge -> (NodeGraph) directedEdge.getFromNode()).collect(Collectors.toList());

    DirectedEdge lastEdge = directedEdgesSequence.get(directedEdgesSequence.size() - 1);
    nodesSequence.add((NodeGraph) lastEdge.getToNode());
  }

  /**
   * Populates the edgesSequence list with edges from the directed edges sequence.
   */
  private void edgeSequence() {

    for (DirectedEdge directedEdge : directedEdgesSequence) {
      edgesSequence.add((EdgeGraph) directedEdge.getEdge());
    }
  }

  /**
   * Computes the sequence of dual nodes for the route using the provided dual network graph.
   *
   * @param dualNetwork The graph representing the dual network.
   */
  public void dualNodesSequence(Graph dualNetwork) {

    directedEdgesSequence.forEach(directedEdge -> {
      EdgeGraph edge = (EdgeGraph) directedEdge.getEdge();
      edgesSequence.add(edge);
      NodeGraph dualNode = dualNetwork.findNode(edge.getDualNode());
      dualNodesSequence.add(dualNode);
    });
  }

  /**
   * Obtains the geometry of the route line by extracting coordinates from the edge sequences.
   */
  private void obtainRouteLineGeometry() {

    Coordinate lastCoordinate = originNode.getCoordinate();
    List<Coordinate> allCoords = new ArrayList<>();

    for (EdgeGraph edge : edgesSequence) {

      LineString geometry = edge.getLine();
      Coordinate[] coords = geometry.getCoordinates();

      // Create a mutable list of coordinates
      List<Coordinate> coordsCollection = new ArrayList<>(Arrays.asList(coords));

      // Check if the direction matches; reverse if necessary
      if (!lastCoordinate.equals(coordsCollection.get(0))) {
        Collections.reverse(coordsCollection);
      }

      // Add the processed coordinates
      allCoords.addAll(coordsCollection);

      // Update the last coordinate
      lastCoordinate = allCoords.get(allCoords.size() - 1);
    }
    // Create a LineString from the combined coordinates
    this.lineString = FACTORY.createLineString(allCoords.toArray(new Coordinate[0]));
  }

  /**
   * Returns the length of the route.
   *
   * @return The length of the route.
   */
  public double getLength() {
    return length;
  }

  /**
   * Returns the LineString representing the geometry of the route.
   *
   * @return A LineString representing the route's geometry.
   */
  public LineString getLineString() {
    return lineString;
  }

  /**
   * Sets the set of visited locations for the route.
   *
   * @param visitedLocations A set of NodeGraph objects representing the visited locations.
   */
  public void setVisitedLocations(Set<NodeGraph> visitedLocations) {
    this.visitedLocations = visitedLocations;
  }

  /**
   * Returns the set of visited locations in the route.
   *
   * @return A set of NodeGraph objects representing the visited locations.
   */
  public Set<NodeGraph> getVisitedLocations() {
    return visitedLocations;
  }
}
