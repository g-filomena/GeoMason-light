/*
 * Copyright (c) 2023 Gabriele Filomena
 * University of Liverpool, UK
 * 
 * This program is free software: it can redistributed and/or modified
 * under the terms of the GNU General Public License 3.0 as published by
 * the Free Software Foundation.
 *
 * See the file "LICENSE" for more information
 */
package sim.graph;

import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.Map;
import java.util.stream.Collectors;

import org.javatuples.Pair;
import org.locationtech.jts.geom.Coordinate;
import org.locationtech.jts.geom.CoordinateArrays;
import org.locationtech.jts.geom.Geometry;
import org.locationtech.jts.geom.LineString;
import org.locationtech.jts.planargraph.DirectedEdge;
import org.locationtech.jts.planargraph.PlanarGraph;

import sim.field.geo.VectorLayer;
import sim.util.geo.MasonGeometry;
import sim.util.geo.Utilities;

/**
 * A planar graph that extends the {@code PlanarGraph} (JTS) class. Its basic components are {@code NodeGraph} and 
 * {@code EdgeGraph}.
 */
public class Graph extends PlanarGraph {
	public ArrayList<EdgeGraph> edgesGraph = new ArrayList<>();
	public ArrayList<NodeGraph> nodesGraph = new ArrayList<>();
	public LinkedHashMap<NodeGraph, Double> centralityMap = new LinkedHashMap<>();
	public HashMap<Integer, NodeGraph> nodesMap = new HashMap<>();
	public VectorLayer junctions = new VectorLayer();
	public HashMap<Pair<NodeGraph, NodeGraph>, EdgeGraph> adjacencyMatrix = new HashMap<>();
	public HashMap<Pair<NodeGraph, NodeGraph>, DirectedEdge> adjacencyMatrixDirected = new HashMap<>();
	public Map<NodeGraph, Double> salientNodes;
	
    /**
     * Constructs a new empty Graph.
     */
	public Graph() {
		super();
	}

	/**
	 * Populates the graph from street junctions and segments. This method adds LineStrings
	 * from the provided street segments to the graph and associates them with nodes from
	 * the street junctions. It also sets the junctions field to the provided street junctions.
	 *
	 * @param streetJunctions The VectorLayer containing street junction geometries.
	 * @param streetSegments  The VectorLayer containing street segment geometries.
	 */
	public void fromStreetJunctionsSegments(VectorLayer streetJunctions, VectorLayer streetSegments) {
		ArrayList<MasonGeometry> geometries = streetSegments.geometriesList;
		for (final MasonGeometry masonGeometry : geometries)
			if (masonGeometry.geometry instanceof LineString)
				addLineString(masonGeometry);
		junctions = streetJunctions;
	}

	/**
	 * Adds a LineString as an EdgeGraph to the graph, creating corresponding nodes if needed.
	 * This method processes the given LineString, removes repeated points, and creates
	 * nodes for the LineString's start and end coordinates. It then constructs an EdgeGraph
	 * representing the LineString and adds it to the graph.
	 *
	 * @param wrappedLine The MasonGeometry object containing the LineString to add.
	 */
	private void addLineString(MasonGeometry wrappedLine) {
		final LineString line = (LineString) wrappedLine.geometry;
		if (line.isEmpty())
			return;

		final Coordinate[] coords = CoordinateArrays.removeRepeatedPoints(line.getCoordinates());
		if (coords.length < 2)
			return;

		final Coordinate fromCoord = coords[0];
		final Coordinate toCoord = coords[coords.length - 1];
		final NodeGraph fromNode = this.getNode(fromCoord);
		final NodeGraph toNode = this.getNode(toCoord);

		nodesMap.put(fromNode.getID(), fromNode);
		nodesMap.put(toNode.getID(), toNode);
		nodesGraph.add(fromNode);
		nodesGraph.add(toNode);

		final EdgeGraph edge = new EdgeGraph(line);
		final DirectedEdge directedEdge0 = new DirectedEdge(fromNode, toNode, coords[1], true);
		final DirectedEdge directedEdge1 = new DirectedEdge(toNode, fromNode, coords[coords.length - 2], false);

		edge.setDirectedEdges(directedEdge0, directedEdge1);
		edge.setAttributes(wrappedLine.getAttributes());
		edge.setNodes(fromNode, toNode);
		edge.masonGeometry = wrappedLine;
		add(edge);
		edgesGraph.add(edge);
	}

	/**
	 * Searches for and retrieves a NodeGraph object in the graph based on the given
	 * coordinate. This method uses the provided coordinate to find a node within the
	 * graph, returning the corresponding NodeGraph object if found, or null if not found.
	 *
	 * @param coordinate The coordinate used to search for the node.
	 * @return The NodeGraph object found at the specified coordinate, or null if not found.
	 */
	@Override
	public NodeGraph findNode(Coordinate coordinate) {
		return (NodeGraph) nodeMap.find(coordinate);
	}

	/**
	 * Retrieves or creates a NodeGraph object based on the given coordinate. If a
	 * node with the specified coordinate already exists in the graph, it is
	 * retrieved. If not, a new NodeGraph object is created, added to the graph, and
	 * returned.
	 *
	 * @param coordinate The coordinate used to find or create the node.
	 * @return The NodeGraph object corresponding to the specified coordinate.
	 */
	public NodeGraph getNode(Coordinate coordinate) {

		NodeGraph node = findNode(coordinate);
		if (node == null) {
			node = new NodeGraph(coordinate);
			// ensure node is only added once to graph
			add(node);
		}
		return node;
	}
	
	/**
	 * Generates structures for the graph. This method generates essential structures, 
	 * including the node mapping (nodesMap) and adjacency matrix (adjacencyMatrix) 
	 * that represent relationships between nodes and edges.
	 */
	public void generateGraphStructures() {
		generateNodesMap();
		generateAdjacencyMatrix();
	}
	
	/**
	 * Finds and stores the salient nodes in the graph based on a given percentile.
	 * This method calculates the salient nodes by selecting nodes with centrality values
	 * above a specified percentile and stores them in the 'salientNodes' map.
	 *
	 * @param salientNodesPercentile The percentile threshold for determining salient nodes.
	 */
	public void setGraphSalientNodes(double salientNodesPercentile) {
		salientNodes = graphSalientNodes(salientNodesPercentile);
	}

	/**
	 * Generates a map of nodes using their IDs as keys. This method creates a
	 * mapping of node IDs to NodeGraph objects for quick and efficient lookup.
	 */
	private void generateNodesMap() {
		for (final NodeGraph node : nodesGraph)
			nodesMap.put(node.getID(), node);
	}

	/**
	 * Generates the centrality map for nodes in the graph. This method computes the
	 * centrality values for each node in the graph and stores them in a LinkedHashMap.
	 * It also rescales the centrality values to the range [0, 1].
	 */
	public void generateCentralityMap() {

		final LinkedHashMap<NodeGraph, Double> centralityMap = new LinkedHashMap<>();
		for (final NodeGraph node : nodesGraph)
			centralityMap.put(node, node.centrality);
		this.centralityMap = (LinkedHashMap<NodeGraph, Double>) Utilities.sortByValue(centralityMap, false);

		// rescale
		for (final NodeGraph node : nodesGraph) {
			final double rescaled = (node.centrality - Collections.min(centralityMap.values()))
					/ (Collections.max(centralityMap.values()) - Collections.min(centralityMap.values()));
			node.centrality_sc = rescaled;
		}
	}
	
	/**
	 * Generates the adjacency matrix and directed adjacency matrix for the graph.
	 * This method populates the adjacency matrix with edges and their corresponding
	 * directed edges, considering both the original and opposite directed edges.
	 * The adjacency matrix stores relationships between nodes and edges.
	 */
	private void generateAdjacencyMatrix() {

		// Populate the adjacency list with edges from your graph
		for (EdgeGraph edgeGraph : edgesGraph) {
			DirectedEdge directedEdge = edgeGraph.getDirEdge(0);
			DirectedEdge oppositeDirectedEdge = edgeGraph.getDirEdge(1);

			// Add the original directed edge and its nodes to the adjacency matrix
			NodeGraph fromNode = (NodeGraph) directedEdge.getFromNode();
			NodeGraph toNode = (NodeGraph) directedEdge.getToNode();
			Pair<NodeGraph, NodeGraph> nodes = new Pair<>(fromNode, toNode);
			adjacencyMatrixDirected.put(nodes, directedEdge);
			adjacencyMatrix.put(nodes, edgeGraph);

			// Add the opposite directed edge and its nodes to the adjacency matrix
			fromNode = (NodeGraph) oppositeDirectedEdge.getFromNode();
			toNode = (NodeGraph) oppositeDirectedEdge.getToNode();
			nodes = new Pair<>(fromNode, toNode);
			adjacencyMatrixDirected.put(nodes, oppositeDirectedEdge);
			adjacencyMatrix.put(nodes, edgeGraph);
		}
	}

	/**
	 * Retrieves the list of edges in the graph.
	 *
	 * @return An ArrayList of EdgeGraph objects representing the edges in the graph.
	 */
	@Override
	public ArrayList<NodeGraph> getNodes() {
		return nodesGraph;
	}

	/**
	 * Retrieves the list of edges in the graph.
	 *
	 * @return An ArrayList of EdgeGraph objects representing the edges in the graph.
	 */
	@Override
	public ArrayList<EdgeGraph> getEdges() {
		return edgesGraph;
	}

	/**
	 * Retrieves a list of graph's nodes contained within the specified Geometry.
	 *
	 * @param geometry The Geometry object used for containment checks.
	 * @return An ArrayList of NodeGraph objects representing nodes that are
	 *         contained within the specified Geometry.
	 */
	public ArrayList<NodeGraph> getContainedNodes(Geometry geometry) {
		final ArrayList<NodeGraph> containedNodes = new ArrayList<>();
		final Collection<NodeGraph> nodes = nodesMap.values();

		for (final NodeGraph node : nodes) {
			final Geometry geoNode = node.masonGeometry.geometry;
			if (geometry.contains(geoNode))
				containedNodes.add(node);
		}
		return containedNodes;
	}

	/**
	 * Retrieves a list of graph's edges contained within the specified Geometry.
	 *
	 * @param geometry The Geometry object used for containment checks.
	 * @return An ArrayList of EdgeGraph objects representing edges that are
	 *         contained within the specified Geometry.
	 */
	public ArrayList<EdgeGraph> getContainedEdges(Geometry geometry) {
		final ArrayList<EdgeGraph> containedEdges = new ArrayList<>();
		final ArrayList<EdgeGraph> edges = edgesGraph;

		for (final EdgeGraph edge : edges) {
			final Geometry edgeGeometry = edge.masonGeometry.geometry;
			if (geometry.contains(edgeGeometry))
				containedEdges.add(edge);
		}
		return containedEdges;
	}

	/**
	 * Retrieves the edge between two nodes, if it exists in the graph's adjacency matrix.
	 *
	 * @param fromNode The source node of the edge.
	 * @param toNode   The target node of the edge.
	 * @return The EdgeGraph object representing the edge between the specified source
	 *         and target nodes, or null if no such edge exists.
	 */
	public EdgeGraph getEdgeBetween(NodeGraph fromNode, NodeGraph toNode) {
		EdgeGraph edge = null;
		Pair<NodeGraph, NodeGraph> pair = new Pair<>(fromNode, toNode);
		edge = adjacencyMatrix.get(pair);
		return edge;
	}

	/**
	 * Retrieves the directed edge between two nodes, if it exists in the graph's
	 * directed adjacency matrix.
	 *
	 * @param fromNode The source node of the directed edge.
	 * @param toNode   The target node of the directed edge.
	 * @return The DirectedEdge object representing the directed edge between the
	 *         specified source and target nodes, or null if no such edge exists.
	 */
	public DirectedEdge getDirectedEdgeBetween(NodeGraph fromNode, NodeGraph toNode) {
		DirectedEdge edge = null;
		Pair<NodeGraph, NodeGraph> pair = new Pair<>(fromNode, toNode);
		edge = adjacencyMatrixDirected.get(pair);
		return edge;
	}

	/**
	 * Filters a LinkedHashMap of centrality values to include only nodes that are present
	 * in the given list of nodes. This method retains entries in the original map for
	 * nodes that are also present in the filter list and returns the filtered map.
	 *
	 * @param map    The LinkedHashMap containing centrality values associated with nodes.
	 * @param filter The list of nodes to use for filtering the map.
	 * @return A filtered LinkedHashMap containing centrality values for nodes present in
	 *         both the original map and the filter list.
	 */
	public static LinkedHashMap<NodeGraph, Double> filterCentralityMap(LinkedHashMap<NodeGraph, Double> map,
			ArrayList<NodeGraph> filter) {

		final LinkedHashMap<NodeGraph, Double> mapFiltered = new LinkedHashMap<>(map);
		final ArrayList<NodeGraph> result = new ArrayList<>();
		for (final NodeGraph key : mapFiltered.keySet())
			if (filter.contains(key))
				result.add(key);
		mapFiltered.keySet().retainAll(result);
		return mapFiltered;
	}

	/**
	 * Retrieves salient nodes from the graph based on a specified percentile threshold
	 * of centrality values. This method calculates the boundary value for centrality
	 * based on the specified percentile and filters nodes with centrality values exceeding
	 * that threshold.
	 *
	 * @param percentile The percentile threshold for selecting salient nodes
	 *                   (e.g., 0.75 for the top 25% centrality).
	 * @return A map of salient nodes and their centrality values exceeding the specified
	 *         percentile threshold.
	 */
	protected Map<NodeGraph, Double> graphSalientNodes(double percentile) {
		int position;
		position = (int) (centralityMap.size() * percentile);
		final double boundary = new ArrayList<>(centralityMap.values()).get(position);

		final Map<NodeGraph, Double> filteredMap = centralityMap.entrySet().stream()
				.filter(entry -> entry.getValue() >= boundary)
				.collect(Collectors.toMap(entry -> entry.getKey(), entry -> entry.getValue()));

		if (filteredMap.isEmpty() || filteredMap == null)
			return null;
		else
			return filteredMap;
	}

	/**
	 * Retrieves salient nodes within the spatial proximity of two reference nodes
	 * based on a specified percentile threshold. This method calculates the smallest
	 * enclosing circle around the reference nodes and filters nodes within that spatial
	 * area. It then filters nodes based on their centrality values and returns a map
	 * of salient nodes exceeding the specified percentile threshold.
	 *
	 * @param node       The first reference node for spatial calculations.
	 * @param otherNode  The second reference node for spatial calculations.
	 * @param percentile The percentile threshold for selecting salient nodes
	 *                   (e.g., 0.75 for the top 25% centrality).
	 * @return A map of salient nodes and their centrality values within the spatial
	 *         proximity of the reference nodes, exceeding the specified percentile
	 *         threshold.
	 */
	public Map<NodeGraph, Double> salientNodesWithinSpace(NodeGraph node, NodeGraph otherNode, double percentile) {

		ArrayList<NodeGraph> containedNodes = new ArrayList<>();
		final Geometry smallestEnclosingCircle = NodeGraphUtils.enclosingCircleBetweenNodes(node, otherNode);
		containedNodes = this.getContainedNodes(smallestEnclosingCircle);

		if (containedNodes.isEmpty())
			return null;
		LinkedHashMap<NodeGraph, Double> spatialfilteredMap = new LinkedHashMap<>();
		spatialfilteredMap = filterCentralityMap(centralityMap, containedNodes);
		if (spatialfilteredMap.isEmpty() || spatialfilteredMap == null)
			return null;

		final int position = (int) (spatialfilteredMap.size() * percentile);
		final double boundary = new ArrayList<>(spatialfilteredMap.values()).get(position);
		final Map<NodeGraph, Double> valueFilteredMap = spatialfilteredMap.entrySet().stream()
				.filter(entry -> entry.getValue() >= boundary)
				.collect(Collectors.toMap(entry -> entry.getKey(), entry -> entry.getValue()));

		if (valueFilteredMap.isEmpty() || valueFilteredMap == null)
			return null;
		else
			return valueFilteredMap;
	}

	/**
	 * Retrieves edges within the spatial proximity of two nodes. This method calculates
	 * a spatial buffer based on the distance between two nodes and retrieves edges
	 * that fall within the combined spatial buffer of the two nodes.
	 *
	 * @param node       The first node used for spatial reference.
	 * @param otherNode  The second node used for spatial reference.
	 * @return A list of edges that are located within the combined spatial buffer of
	 *         the two nodes.
	 */
	public ArrayList<EdgeGraph> edgesInNodesSpace(NodeGraph node, NodeGraph otherNode) {

		Double radius = NodeGraphUtils.nodesDistance(node, otherNode) * 1.50;
		if (radius < 500)
			radius = 500.0;
		final Geometry bufferOrigin = node.masonGeometry.geometry.buffer(radius);
		final Geometry bufferDestination = otherNode.masonGeometry.geometry.buffer(radius);
		final Geometry convexHull = bufferOrigin.union(bufferDestination).convexHull();

		final ArrayList<EdgeGraph> containedEdges = this.getContainedEdges(convexHull);
		return containedEdges;
	}

	/**
	 * Filters a list of nodes to exclude nodes belonging to a specific region.
	 *
	 * @param nodes    The list of nodes to filter.
	 * @param regionID The ID of the region to exclude nodes from.
	 * @return A new list containing nodes that do not belong to the specified region.
	 */
	public ArrayList<NodeGraph> nodesInRegion(ArrayList<NodeGraph> nodes, int regionID) {
		final ArrayList<NodeGraph> newNodes = new ArrayList<>(nodes);
		for (final NodeGraph node : nodes)
			if (node.regionID == regionID)
				newNodes.remove(node);
		return newNodes;
	}
}
