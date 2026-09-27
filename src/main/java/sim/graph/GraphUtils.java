/*
 * Copyright (c) 2023 Gabriele Filomena University of Liverpool, UK
 *
 * This program is free software: it can redistributed and/or modified under the terms of the GNU
 * General Public License 3.0 as published by the Free Software Foundation.
 *
 * See the file "LICENSE" for more information
 */
package sim.graph;

import java.util.ArrayList;
import java.util.Comparator;
import java.util.List;
import java.util.Map;
import java.util.Objects;
import java.util.Set;
import java.util.stream.Collectors;
import java.util.stream.StreamSupport;
import org.locationtech.jts.algorithm.ConvexHull;
import org.locationtech.jts.geom.Coordinate;
import org.locationtech.jts.geom.Geometry;
import org.locationtech.jts.geom.GeometryFactory;
import org.locationtech.jts.geom.LineString;
import org.locationtech.jts.geom.Point;
import org.locationtech.jts.geom.Polygon;
import sim.util.geo.GeometryUtilities;

/**
 * The `GraphUtils` class provides utility methods for working {@link NodeGraph} objects.
 */
public class GraphUtils {

  private static final GeometryFactory GEOMETRY_FACTORY = new GeometryFactory();

  /**
   * Calculates the Euclidean distance between two nodes in a graph.
   *
   * @param node The first node.
   * @param otherNode The second node.
   * @return The Euclidean distance between the two nodes' coordinates.
   */
  public static double nodesDistance(NodeGraph node, NodeGraph otherNode) {

    // Computed, not cached. The static cache this replaces kept every pair ever measured for the
    // life of the JVM - an unbounded leak in a long run - and a square root is cheaper than the two
    // Pair allocations and hash lookups it took to consult it.
    return GeometryUtilities.euclideanDistance(node.getCoordinate(), otherNode.getCoordinate());
  }

  /**
   * Calculates the minimum enclosing circle (smallest circle that completely encloses a collection
   * of nodes).
   *
   * @param nodes The collection of nodes for which to calculate the enclosing circle.
   * @return A geometry representing the minimum enclosing circle.
   */
  public static Geometry smallestEnclosingGeometryBetweenNodes(List<NodeGraph> nodes) {

    if (nodes.size() == 1) {
      return nodes.get(0).masonGeometry.getGeometry().buffer(50);
    }
    if (nodes.size() == 2) {
      return enclosingCircleBetweenTwoNodes(nodes.get(0), nodes.get(1));
    }
    return convexHullFromNodes(nodes);
  }

  /**
   * Calculates the convex hull of a list of nodes.
   *
   * <p>Nodes lying on one line have a hull with no area, which contains none of them; for those
   * this returns the circle through the two outermost nodes, as for two nodes, and for nodes all
   * at one point the same 50-unit buffer as for a single node. The caller's list is not reordered.
   *
   * @param nodes The list of nodes to compute the convex hull from.
   * @return A geometry enclosing all the nodes.
   */
  private static Geometry convexHullFromNodes(List<NodeGraph> nodes) {

    Coordinate[] coordinates = nodes.stream().map(NodeGraph::getCoordinate)
        .toArray(Coordinate[]::new);
    Geometry hull = new ConvexHull(coordinates, GEOMETRY_FACTORY).getConvexHull();
    if (hull instanceof Polygon) {
      return hull;
    }
    if (hull instanceof LineString) {
      // For collinear input the hull is the segment between the two outermost points.
      LineString segment = (LineString) hull;
      return enclosingCircle(segment.getCoordinateN(0),
          segment.getCoordinateN(segment.getNumPoints() - 1));
    }
    return hull.buffer(50);
  }

  /**
   * Calculates the smallest enclosing circle for two nodes in a graph.
   *
   * @param node The first node.
   * @param otherNode The second node.
   * @return The smallest enclosing circle as a geometry.
   */
  protected static Geometry enclosingCircleBetweenTwoNodes(NodeGraph node, NodeGraph otherNode) {
    return enclosingCircle(node.getCoordinate(), otherNode.getCoordinate());
  }

  /** The circle whose diameter is the segment between the two coordinates. */
  private static Geometry enclosingCircle(Coordinate coordinate, Coordinate otherCoordinate) {
    final LineString line =
        GEOMETRY_FACTORY.createLineString(new Coordinate[] {coordinate, otherCoordinate});
    final Point centroid = line.getCentroid();
    return centroid.buffer(line.getLength() / 2);
  }

  /**
   * Finds the closest node to a target coordinate from a collection of nodes.
   *
   * @param targetCoordinates The target coordinate to which to find the closest node.
   * @param nodes The collection of nodes to search for the closest node.
   * @return The closest node to the target coordinate.
   */
  public static NodeGraph findClosestNode(Coordinate targetCoordinates, Iterable<NodeGraph> nodes) {
    return StreamSupport.stream(nodes.spliterator(), true) // Convert Iterable to parallel stream
        .min(Comparator.comparingDouble(node -> targetCoordinates.distance(node.getCoordinate())))
        .orElse(null); // found
  }

  /**
   * It returns a LineString between two given nodes.
   *
   * @param node a node;
   * @param otherNode an other node;
   */
  public static LineString LineStringBetweenNodes(NodeGraph node, NodeGraph otherNode) {
    final Coordinate[] coords = {node.getCoordinate(), otherNode.getCoordinate()};
    final LineString line = new GeometryFactory().createLineString(coords);
    return line;
  }

  /**
   * Extracts all nodes from a given set of edges.
   *
   * @param edges the set of edges from which to extract nodes
   * @return a set containing all unique nodes from the given edges
   */
  public static Set<NodeGraph> nodesFromEdges(Set<EdgeGraph> edges) {
    // Insertion-ordered: Islands seeds its search from these in turn, and NodeGraph inherits
    // identity hashing, which differs between JVM builds. Collectors.toSet() promises no order.
    return edges.stream().flatMap(edge -> edge.getNodes().stream())
        .collect(Collectors.toCollection(java.util.LinkedHashSet::new));
  }

  /**
   * Extracts all edges from a given set of nodes.
   *
   * @param nodes the set of nodes from which to extract edges
   * @return a set containing all unique edges from the given nodes
   */
  public static Set<EdgeGraph> edgesFromNodes(Set<NodeGraph> nodes) {
    return nodes.stream().flatMap(node -> node.getEdges().stream())
        .collect(Collectors.toCollection(java.util.LinkedHashSet::new));
  }

  /**
   * Retrieves the IDs of a given list of nodes.
   *
   * @param nodes the list of nodes from which to extract IDs
   * @return a list of IDs corresponding to the given nodes
   */
  public static List<Integer> getNodeIDs(List<NodeGraph> nodes) {
    return nodes.stream().map(NodeGraph::getID).collect(Collectors.toList());
  }

  /**
   * Retrieves the IDs of a given list of nodes.
   *
   * @param nodes the list of nodes from which to extract IDs
   * @return a list of IDs corresponding to the given nodes
   */
  public static List<Integer> getNodeIDs(Set<NodeGraph> nodes) {
    return nodes.stream().map(NodeGraph::getID).collect(Collectors.toList());
  }

  /**
   * Retrieves the IDs of a given list of edges.
   *
   * @param edges the list of edges from which to extract IDs
   * @return a list of IDs corresponding to the given edges
   */
  public static List<Integer> getEdgeIDs(List<EdgeGraph> edges) {
    return edges.stream().map(EdgeGraph::getID).collect(Collectors.toList());
  }

  /**
   * Retrieves the IDs of a given list of edges.
   *
   * @param edges the list of edges from which to extract IDs
   * @return a list of IDs corresponding to the given edges
   */
  public static List<Integer> getEdgeIDs(Set<EdgeGraph> edges) {
    return edges.stream().map(EdgeGraph::getID).collect(Collectors.toList());
  }

  /**
   * Retrieves nodes from a list of node IDs using a provided map.
   *
   * @param nodeIDs A list of node IDs.
   * @param map The map to retrieve nodes from.
   * @return A list of nodes corresponding to the given node IDs.
   */
  public static List<NodeGraph> getNodesFromNodeIDs(List<Integer> nodeIDs,
      Map<Integer, NodeGraph> map) {
    return nodeIDs.stream().map(map::get) // Retrieve NodeGraph from the map
        .filter(Objects::nonNull) // Ignore null values
        .collect(Collectors.toList());
  }

  /**
   * Retrieves edges from a list of edge IDs using a provided map.
   *
   * @param edgeIDs A list of edge IDs.
   * @param map The map to retrieve edges from.
   * @return A list of edges corresponding to the given edge IDs.
   */
  public static List<EdgeGraph> getEdgesFromEdgeIDs(List<Integer> edgeIDs,
      Map<Integer, EdgeGraph> map) {
    return edgeIDs.stream().map(map::get) // Retrieve EdgeGraph from the map
        .filter(Objects::nonNull) // Ignore null values
        .collect(Collectors.toList());
  }

  /**
   * Retrieves nodes from a set of node IDs using a provided map.
   *
   * @param nodeIDs A set of node IDs.
   * @param map The map to retrieve nodes from.
   * @return A list of nodes corresponding to the given node IDs.
   */
  public static List<NodeGraph> getNodesFromNodeIDs(Set<Integer> nodeIDs,
      Map<Integer, NodeGraph> map) {
    return nodeIDs.stream().map(map::get) // Retrieve NodeGraph from the map
        .filter(Objects::nonNull) // Ignore null values
        .collect(Collectors.toList());
  }

  /**
   * Retrieves edges from a set of edge IDs using a provided map.
   *
   * @param edgeIDs A set of edge IDs.
   * @param map The map to retrieve edges from.
   * @return A list of edges corresponding to the given edge IDs.
   */
  public static List<EdgeGraph> getEdgesFromEdgeIDs(Set<Integer> edgeIDs,
      Map<Integer, EdgeGraph> map) {
    return edgeIDs.stream().map(map::get) // Retrieve EdgeGraph from the map
        .filter(Objects::nonNull) // Ignore null values
        .collect(Collectors.toList());
  }
}
