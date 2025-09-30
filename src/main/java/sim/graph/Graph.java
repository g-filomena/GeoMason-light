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
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;
import org.javatuples.Pair;
import org.locationtech.jts.geom.Coordinate;
import org.locationtech.jts.geom.CoordinateArrays;
import org.locationtech.jts.geom.Geometry;
import org.locationtech.jts.geom.GeometryFactory;
import org.locationtech.jts.geom.LineString;
import org.locationtech.jts.geom.Polygon;
import org.locationtech.jts.index.strtree.STRtree;
import org.locationtech.jts.planargraph.DirectedEdge;
import org.locationtech.jts.planargraph.PlanarGraph;
import sim.field.geo.VectorLayer;
import sim.util.geo.AttributeValue;
import sim.util.geo.MasonGeometry;
import sim.util.geo.Utilities;

/**
 * Represents a planar graph that extends the JTS (Java Topology Suite) PlanarGraph class. Its basic
 * components are {@link NodeGraph} and {@link EdgeGraph}. This class provides functionalities to
 * construct and manipulate a graph from street junctions and segments.
 */
public class Graph extends PlanarGraph {

  protected List<NodeGraph> nodesGraph = new ArrayList<>();
  protected List<EdgeGraph> edgesGraph = new ArrayList<>();
  protected Map<NodeGraph, Double> centralityMap = new LinkedHashMap<>();
  protected VectorLayer junctions = new VectorLayer();

  protected Map<Pair<NodeGraph, NodeGraph>, EdgeGraph> adjacencyMatrix = new HashMap<>();
  protected Map<Pair<NodeGraph, NodeGraph>, DirectedEdge> adjacencyMatrixDirected = new HashMap<>();
  protected Map<NodeGraph, Double> salientNodes = new HashMap<>();
  protected Set<DirectedEdge> directedEdges = new HashSet<>();
  public Map<String, AttributeValue> attributes = new HashMap<>();

  // Build spatial index once
  STRtree index = new STRtree();

  private static final GeometryFactory GEOMETRY_FACTORY = new GeometryFactory();

  /**
   * Constructs a new empty Graph.
   */
  public Graph() {
    super();
  }

  /**
   * Populates the graph with nodes and edges based on street junctions and segments. Adds
   * LineStrings from the street segments to the graph and creates nodes at the junctions. Also sets
   * the 'junctions' field with the provided street junctions.
   *
   * @param streetJunctions The VectorLayer containing street junction geometries.
   * @param streetSegments The VectorLayer containing street segment geometries.
   */
  public void fromStreetJunctionsSegments(VectorLayer streetJunctions, VectorLayer streetSegments) {

    List<MasonGeometry> geometries = streetSegments.getGeometries();
    geometries.stream().filter(masonGeometry -> masonGeometry.geometry instanceof LineString)
        .forEach(this::addLineString);

    nodesGraph.forEach(NodeGraph::setNeighbouringComponents);
    generateAdjacencyMatrix();

    for (NodeGraph node : getNodes()) {
      index.insert(node.masonGeometry.getGeometry().getEnvelopeInternal(), node);
    }
  }

  /**
   * Populates the graph with nodes and edges based on street junctions and segments. Adds
   * LineStrings from the street segments to the graph and creates nodes at the junctions. Also sets
   * the 'junctions' field with the provided street junctions.
   *
   * @param streetJunctions The VectorLayer containing street junction geometries.
   * @param streetSegments The VectorLayer containing street segment geometries.
   */
  private void addLineString(MasonGeometry wrappedLine) {
    final LineString line = (LineString) wrappedLine.geometry;
    if (line.isEmpty()) {
      return;
    }

    Coordinate[] coords = CoordinateArrays.removeRepeatedPoints(line.getCoordinates());
    if (coords.length < 2) {
      return;
    }

    Coordinate fromCoord = coords[0];
    Coordinate toCoord = coords[coords.length - 1];
    NodeGraph fromNode = getNode(fromCoord);
    NodeGraph toNode = getNode(toCoord);

    EdgeGraph edge = new EdgeGraph(line);
    DirectedEdge directedEdge0 = new DirectedEdge(fromNode, toNode, coords[1], true);
    DirectedEdge directedEdge1 =
        new DirectedEdge(toNode, fromNode, coords[coords.length - 2], false);

    edge.setDirectedEdges(directedEdge0, directedEdge1);
    edge.setAttributes(wrappedLine.getAttributes());
    edge.setNodes(fromNode, toNode);
    edge.masonGeometry = wrappedLine;
    addEdge(edge);
  }

  /**
   * Adds an EdgeGraph to the graph. This method includes the EdgeGraph in the list of graph edges
   * and adds it to the PlanarGraph structure. It ensures the edge is part of the graph's edge
   * collection and its planar graph representation.
   *
   * @param edge The EdgeGraph to be added to the graph.
   */
  protected void addEdge(EdgeGraph edge) {
    edgesGraph.add(edge);
    add(edge);
  }

  /**
   * Searches for and retrieves a NodeGraph object in the graph based on the given coordinate. This
   * method uses the provided coordinate to find a node within the graph, returning the
   * corresponding NodeGraph object if found, or null if not found.
   *
   * @param coordinate The coordinate used to search for the node.
   * @return The NodeGraph object found at the specified coordinate, or null if not found.
   */
  @Override
  public NodeGraph findNode(Coordinate coordinate) {
    return (NodeGraph) nodeMap.find(coordinate);
  }

  /**
   * Retrieves or creates a NodeGraph object based on the given coordinate. If a node with the
   * specified coordinate already exists in the graph, it is retrieved. If not, a new NodeGraph
   * object is created, added to the graph, and returned.
   *
   * @param coordinate The coordinate used to find or create the node.
   * @return The NodeGraph object corresponding to the specified coordinate.
   */
  public NodeGraph getNode(Coordinate coordinate) {

    NodeGraph node = findNode(coordinate);
    if (node == null) {
      node = new NodeGraph(coordinate);

      // ensure node is only added once to graph
      node.masonGeometry = new MasonGeometry(GEOMETRY_FACTORY.createPoint(coordinate));
      nodesGraph.add(node);
      add(node);
      junctions.addGeometry(node.masonGeometry);
    }
    return node;
  }

  /**
   * Generates the centrality map for nodes in the graph. This method computes the centrality values
   * for each node in the graph and stores them in a LinkedHashMap.
   */
  public void generateCentralityMap() {

    Map<NodeGraph, Double> centralityMap = new LinkedHashMap<>();
    for (NodeGraph node : nodesGraph) {
      centralityMap.put(node, node.centrality);
    }
    this.centralityMap = Utilities.sortByValue(centralityMap, false);
  }

  /**
   * Generates the adjacency matrix and directed adjacency matrix for the graph. This method
   * populates the adjacency matrix with edges and their corresponding directed edges, considering
   * both the original and opposite directed edges. The adjacency matrix stores relationships
   * between nodes and edges.
   */
  protected void generateAdjacencyMatrix() {
    for (EdgeGraph edgeGraph : edgesGraph) {
      addEdgesToAdjacencyMatrix(edgeGraph, 0);
      addEdgesToAdjacencyMatrix(edgeGraph, 1);
    }
  }

  private void addEdgesToAdjacencyMatrix(EdgeGraph edgeGraph, int index) {
    DirectedEdge directedEdge = edgeGraph.getDirEdge(index);
    NodeGraph fromNode = (NodeGraph) directedEdge.getFromNode();
    NodeGraph toNode = (NodeGraph) directedEdge.getToNode();
    Pair<NodeGraph, NodeGraph> nodes = new Pair<>(fromNode, toNode);

    adjacencyMatrixDirected.put(nodes, directedEdge);
    adjacencyMatrix.put(nodes, edgeGraph);
  }

  /**
   * Retrieves the list of edges in the graph.
   *
   * @return A List of EdgeGraph objects representing the edges in the graph.
   */
  @Override
  public List<NodeGraph> getNodes() {
    return nodesGraph;
  }

  /**
   * Retrieves the IDs of all nodes in the graph.
   *
   * @return a list of IDs corresponding to the nodes in the graph
   */
  public List<Integer> getNodeIDs() {

    List<Integer> nodeIDs = new ArrayList<>();
    for (NodeGraph node : nodesGraph) {
      nodeIDs.add(node.getID());
    }
    return nodeIDs;
  }

  /**
   * Retrieves the IDs of all edges in the graph.
   *
   * @return a list of IDs corresponding to the edges in the graph
   */
  public List<Integer> getEdgeIDs() {

    List<Integer> edgeIDs = new ArrayList<>();
    for (EdgeGraph edge : edgesGraph) {
      edgeIDs.add(edge.getID());
    }
    return edgeIDs;
  }

  /**
   * Retrieves the list of edges in the graph.
   *
   * @return A List of EdgeGraph objects representing the edges in the graph.
   */
  @Override
  public List<EdgeGraph> getEdges() {
    return edgesGraph;
  }

  /**
   * Retrieves a List of graph's nodes contained within the specified Geometry.
   *
   * @param geometry The Geometry object used for containment checks.
   * @return A List of NodeGraph objects representing nodes that are contained within the specified
   *         Geometry.
   */
  public List<NodeGraph> getContainedNodes(Geometry geometry) {

    return nodesGraph.parallelStream()
        .filter(node -> geometry.contains(node.masonGeometry.geometry))
        .collect(Collectors.toList());
  }

  /**
   * Retrieves a List of graph's edges contained within the specified Geometry.
   *
   * @param geometry The Geometry object used for containment checks.
   * @return A List of EdgeGraph objects representing edges that are contained within the specified
   *         Geometry.
   */
  public List<EdgeGraph> getContainedEdges(Geometry geometry) {
    return edgesGraph.parallelStream()
        .filter(edge -> geometry.contains(edge.masonGeometry.geometry))
        .collect(Collectors.toList());
  }

  /**
   * Retrieves the edge between two nodes, if it exists in the graph's adjacency matrix.
   *
   * @param fromNode The source node of the edge.
   * @param toNode The target node of the edge.
   * @return The EdgeGraph object representing the edge between the specified source and target
   *         nodes, or null if no such edge exists.
   */
  public EdgeGraph getEdgeBetween(NodeGraph fromNode, NodeGraph toNode) {
    return adjacencyMatrix.get(new Pair<>(fromNode, toNode));
  }

  /**
   * Retrieves the directed edge between two nodes, if it exists in the graph's directed adjacency
   * matrix.
   *
   * @param fromNode The source node of the directed edge.
   * @param toNode The target node of the directed edge.
   * @return The DirectedEdge object representing the directed edge between the specified source and
   *         target nodes, or null if no such edge exists.
   */
  public DirectedEdge getDirectedEdgeBetween(NodeGraph fromNode, NodeGraph toNode) {
    return adjacencyMatrixDirected.get(new Pair<>(fromNode, toNode));
  }

  /**
   * Filters a LinkedHashMap of centrality values to include only nodes that are present in the
   * given list of nodes. This method retains entries in the original map for nodes that are also
   * present in the filter list and returns the filtered map.
   *
   * @param map The LinkedHashMap containing centrality values associated with nodes.
   * @param filter The list of nodes to use for filtering the map.
   * @return A filtered LinkedHashMap containing centrality values for nodes present in both the
   *         original map and the filter list.
   */
  public static Map<NodeGraph, Double> filterCentralityMap(Map<NodeGraph, Double> map,
      List<NodeGraph> filter) {

    return map.entrySet().stream().filter(entry -> filter.contains(entry.getKey()))
        .collect(Collectors.toMap(Map.Entry::getKey, Map.Entry::getValue, (e1, e2) -> e1,
            LinkedHashMap::new));
  }

  /**
   * Retrieves salient nodes from the graph based on a specified percentile threshold of centrality
   * values. This method calculates the boundary value for centrality based on the specified
   * percentile and filters nodes with centrality values exceeding that threshold.
   *
   * @param percentile The percentile threshold for selecting salient nodes (e.g., 0.75 for the top
   *        25% centrality).
   * @return A map of salient nodes and their centrality values exceeding the specified percentile
   *         threshold.
   */
  public Map<NodeGraph, Double> getSalientNodes(double percentile) {

    int position = (int) (centralityMap.size() * percentile);
    ArrayList<Double> values = new ArrayList<>(centralityMap.values());
    values.sort(Double::compareTo); // Sort values to determine the boundary
    double boundary = values.get(position);

    return centralityMap.entrySet().stream().filter(entry -> entry.getValue() >= boundary)
        .collect(Collectors.toMap(Map.Entry::getKey, Map.Entry::getValue));

  }

  /**
   * Retrieves salient nodes within the spatial proximity of two reference nodes based on a
   * specified percentile threshold. This method calculates the smallest enclosing circle around the
   * reference nodes and filters nodes within that spatial area. It then filters nodes based on
   * their centrality values and returns a map of salient nodes exceeding the specified percentile
   * threshold.
   *
   * @param node The first reference node for spatial calculations.
   * @param otherNode The second reference node for spatial calculations.
   * @param percentile The percentile threshold for selecting salient nodes (e.g., 0.75 for the top
   *        25% centrality).
   * @return A map of salient nodes and their centrality values within the spatial proximity of the
   *         reference nodes, exceeding the specified percentile threshold.
   */
  public Map<NodeGraph, Double> getSalientNodesWithinSpace(NodeGraph node, NodeGraph otherNode,
      double percentile) {

    List<NodeGraph> containedNodes = new ArrayList<>();
    Geometry smallestEnclosingCircle = GraphUtils.enclosingCircleBetweenTwoNodes(node, otherNode);
    containedNodes = this.getContainedNodes(smallestEnclosingCircle);

    Map<NodeGraph, Double> spatialfilteredMap = new LinkedHashMap<>();
    if (containedNodes.isEmpty()) {
      return spatialfilteredMap;
    }
    spatialfilteredMap = filterCentralityMap(new LinkedHashMap<>(centralityMap), containedNodes);

    if (spatialfilteredMap.isEmpty()) {
      return spatialfilteredMap;
    }

    int position = (int) (spatialfilteredMap.size() * percentile);
    double boundary = new ArrayList<>(spatialfilteredMap.values()).get(position);
    Map<NodeGraph, Double> valueFilteredMap =
        spatialfilteredMap.entrySet().stream().filter(entry -> entry.getValue() >= boundary)
            .collect(Collectors.toMap(entry -> entry.getKey(), entry -> entry.getValue()));

    return valueFilteredMap;
  }

  /**
   * Returns a List of nodes in the graph whose distance from the given node falls within the
   * specified range and have centrality values above the specified percentile.
   *
   * @param node The reference node.
   * @param lowerLimit The minimum distance from the reference node.
   * @param upperLimit The maximum distance from the reference node.
   * @param percentile The percentile used as a threshold for centrality values.
   * @return A List of nodes that meet the distance and centrality criteria.
   */
  public List<NodeGraph> getSalientNodesBetweenDistanceInterval(NodeGraph node, double lowerLimit,
      double upperLimit, double percentile) {

    List<NodeGraph> containedNodes = new ArrayList<>();
    List<NodeGraph> containedSalientNodes = new ArrayList<>();
    MasonGeometry originGeometry = node.masonGeometry;
    List<MasonGeometry> containedGeometries =
        junctions.featuresBetweenLimits(originGeometry.geometry, lowerLimit, upperLimit);
    for (MasonGeometry masonGeometry : containedGeometries) {
      containedNodes.add(findNode(masonGeometry.geometry.getCoordinate()));
    }

    List<NodeGraph> salientNodes = new ArrayList<>(getSalientNodes(percentile).keySet());

    for (NodeGraph otherNode : containedNodes) {
      if (salientNodes.contains(otherNode)) {
        containedSalientNodes.add(otherNode);
      }
    }
    return containedSalientNodes;
  }

  /**
   * Retrieves edges within the spatial proximity of two nodes. This method calculates a spatial
   * buffer based on the distance between two nodes and retrieves edges that fall within the
   * combined spatial buffer of the two nodes.
   *
   * @param node The first node used for spatial reference.
   * @param otherNode The second node used for spatial reference.
   * @return A List of edges that are located within the combined spatial buffer of the two nodes.
   */
  public List<EdgeGraph> edgesInNodesSpace(NodeGraph node, NodeGraph otherNode) {

    double radius = Math.max(GraphUtils.nodesDistance(node, otherNode) * 1.5, 500.0);
    Geometry bufferOrigin = node.masonGeometry.geometry.buffer(radius);
    Geometry bufferDestination = otherNode.masonGeometry.geometry.buffer(radius);
    Geometry convexHull = bufferOrigin.union(bufferDestination).convexHull();
    return getContainedEdges(convexHull);
  }

  /**
   * Filters out nodes belonging to a specified region from a List of nodes. Returns a new list
   * containing nodes not associated with the specified region ID.
   *
   * @param nodes The original list of nodes to filter.
   * @param regionID The ID of the region to exclude nodes from.
   * @return A List of nodes not belonging to the specified region.
   */
  public List<NodeGraph> nodesInRegion(List<NodeGraph> nodes, int regionID) {
    return nodes.stream().filter(node -> node.regionID != regionID).collect(Collectors.toList());
  }

  /**
   * Finds a node in the graph that matches the given node's geometry coordinates.
   *
   * @param nodeGraph the node whose geometry coordinates will be used for the search
   * @return the node in the graph that matches the given node's geometry coordinates
   */
  public NodeGraph findNode(NodeGraph nodeGraph) {
    return findNode(nodeGraph.getMasonGeometry().geometry.getCoordinate());
  }

  /**
   * Retrieves all {@link NodeGraph} objects whose geometries are contained within the specified
   * polygon.
   *
   * The method first queries the spatial index with the polygon's bounding envelope to retrieve
   * candidate nodes. If no candidates are found, an empty list is returned. Otherwise, a precise
   * containment check is performed against the given polygon to filter only those nodes whose
   * geometry is strictly contained within it.
   *
   * @param polygon the {@link Polygon} defining the search area; must not be null
   * @return a list of {@link NodeGraph} instances located inside the polygon; the list is empty if
   *         no nodes fall within the polygon
   */
  public List<NodeGraph> getNodesWithinPolygon(Polygon polygon) {
    List<NodeGraph> nodesInPolygon = new ArrayList<>();
    @SuppressWarnings("unchecked")
    List<NodeGraph> candidates =
        (List<NodeGraph>) (List<?>) index.query(polygon.getEnvelopeInternal());

    if (candidates.isEmpty()) {
      return nodesInPolygon;
    }

    return candidates.stream().filter(node -> polygon.contains(node.masonGeometry.getGeometry()))
        .collect(Collectors.toList());
  }
}
