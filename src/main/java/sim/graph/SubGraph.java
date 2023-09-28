/* 
 * Copyright 2023 by Gabriele Filomena
 * University of Liverpool, UK
 * The MIT License (MIT)
 *
 */

package sim.graph;

import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

import org.locationtech.jts.geom.Coordinate;
import org.locationtech.jts.geom.CoordinateArrays;
import org.locationtech.jts.geom.LineString;
import org.locationtech.jts.planargraph.DirectedEdge;

import sim.util.geo.Utilities;

/**
 * This class represents sub graphs derived from a graph. It establishes links
 * between the components of the parent graph and the corresponding child
 * components, allowing faster operations and the creation of "regional" or
 * "district" graphs.
 *
 * Agent navigation throughout a SubGraph is straightforward and can be easily
 * retraced to the parent graph.
 */

public class SubGraph extends Graph {
	private final SubGraphNodesMap subGraphNodesMap = new SubGraphNodesMap();
	private final SubGraphEdgesMap subGraphEdgesMap = new SubGraphEdgesMap();
	private ArrayList<Integer> graphBarriers = new ArrayList<>();
	LinkedHashMap<NodeGraph, Double> centralityMap = new LinkedHashMap<>();
	Graph parentGraph = new Graph();

	/**
	 * Constructs a subgraph from a parent graph by copying a specified list of edges.
	 *
	 * @param parentGraph The parent graph from which the edges are copied.
	 * @param edges       The list of edges to be included in the subgraph.
	 */
	public SubGraph(Graph parentGraph, ArrayList<EdgeGraph> edges) {
	    this.parentGraph = parentGraph;
	    for (final EdgeGraph edge : edges)
	        this.addFromOtherGraph(parentGraph, edge);
	    for (final NodeGraph node : this.getNodesList())
	        node.setNeighbouringComponents();
	    generateGraphStructures();
	}

	/**
	 * Constructs an empty subgraph.
	 */
	public SubGraph() {
	}

	/**
	 * Adds an edge and its corresponding nodes to the subgraph from the parent graph.
	 *
	 * @param parentGraph The parent graph from which the edge is copied.
	 * @param parentEdge  The edge to be added to the subgraph.
	 */
	public void addFromOtherGraph(Graph parentGraph, EdgeGraph parentEdge) {
	    final NodeGraph fromNode = parentEdge.fromNode;
	    final NodeGraph toNode = parentEdge.toNode;
	    final Coordinate fromNodeCoord = fromNode.getCoordinate();
	    final Coordinate toNodeCoord = toNode.getCoordinate();

	    final NodeGraph childFromNode = getNode(fromNodeCoord);
	    final NodeGraph childToNode = getNode(toNodeCoord);
	    childFromNode.nodeID = fromNode.getID();
	    childToNode.nodeID = fromNode.getID();

	    final LineString line = parentEdge.getLine();
	    final Coordinate[] coords = CoordinateArrays.removeRepeatedPoints(line.getCoordinates());
	    final EdgeGraph childEdge = new EdgeGraph(line);

	    final DirectedEdge de0 = new DirectedEdge(childFromNode, childToNode, coords[1], true);
	    final DirectedEdge de1 = new DirectedEdge(childToNode, childFromNode, coords[coords.length - 2], false);
	    childEdge.setDirectedEdges(de0, de1);
	    childEdge.setNodes(childFromNode, childToNode);
	    setAttributesChildEdge(childEdge, parentEdge);

	    subGraphNodesMap.add(childFromNode, fromNode);
	    subGraphNodesMap.add(childToNode, toNode);
	    subGraphEdgesMap.add(childEdge, parentEdge);
	    childFromNode.primalEdge = fromNode.primalEdge;
	    childToNode.primalEdge = toNode.primalEdge;
	    add(childEdge);
	    edgesGraph.add(childEdge);
	}

	/**
	 * A mapping utility class for managing relationships between nodes in a subgraph.
	 */
	private class SubGraphNodesMap {
	    public HashMap<NodeGraph, NodeGraph> map = new HashMap<>();

	    /**
	     * Adds a mapping between a child node and its parent node in the subgraph.
	     *
	     * @param node       The child node to be added to the subgraph.
	     * @param parentNode The parent node corresponding to the child node.
	     */
	    private void add(NodeGraph node, NodeGraph parentNode) {
	        map.put(node, parentNode);
	    }

	    /**
	     * Retrieves the parent node corresponding to a given child node in the subgraph.
	     *
	     * @param nodeSubGraph The child node for which the parent node is to be retrieved.
	     * @return The parent node corresponding to the provided child node, or null if not found.
	     */
	    private NodeGraph findParent(NodeGraph nodeSubGraph) {
	        return map.get(nodeSubGraph);
	    }

	    /**
	     * Retrieves the child node corresponding to a given parent node in the subgraph.
	     *
	     * @param nodeGraph The parent node for which the child node is to be retrieved.
	     * @return The child node corresponding to the provided parent node, or null if not found.
	     */
	    private NodeGraph findChild(NodeGraph nodeGraph) {
	        return Utilities.getKeyFromValue(map, nodeGraph);
	    }
	}

	/**
	 * A mapping utility class for managing relationships between edges in a subgraph.
	 */
	private class SubGraphEdgesMap {

	    private final HashMap<EdgeGraph, EdgeGraph> map = new HashMap<>();

	    /**
	     * Adds a mapping between a child edge and its parent edge in the subgraph.
	     *
	     * @param edge        The child edge to be added to the subgraph.
	     * @param parentEdge  The parent edge corresponding to the child edge.
	     */
	    public void add(EdgeGraph edge, EdgeGraph parentEdge) {
	        this.map.put(edge, parentEdge);
	    }

	    /**
	     * Retrieves the parent edge corresponding to a given child edge in the subgraph.
	     *
	     * @param edgeSubGraph The child edge for which the parent edge is to be retrieved.
	     * @return The parent edge corresponding to the provided child edge, or null if not found.
	     */
	    private EdgeGraph findParent(EdgeGraph edgeSubGraph) {
	        return this.map.get(edgeSubGraph);
	    }

	    /**
	     * Retrieves the child edge corresponding to a given parent edge in the subgraph.
	     *
	     * @param edgeSubGraph The parent edge for which the child edge is to be retrieved.
	     * @return The child edge corresponding to the provided parent edge, or null if not found.
	     */
	    private EdgeGraph findChild(EdgeGraph edgeSubGraph) {
	        return Utilities.getKeyFromValue(this.map, edgeSubGraph);
	    }
	}

	/**
	 * Retrieves the parent node corresponding to the provided child node.
	 *
	 * @param childNode The child node for which the parent node is to be retrieved.
	 * @return The parent node corresponding to the provided child node, or null if not found.
	 */
	public NodeGraph getParentNode(NodeGraph childNode) {
	    return this.subGraphNodesMap.findParent(childNode);
	}

	/**
	 * Retrieves a list of all parent nodes within the current subgraph.
	 *
	 * @return A list of all parent nodes present in the current subgraph.
	 */
	public ArrayList<NodeGraph> getParentNodes() {
	    final ArrayList<NodeGraph> parentNodes = new ArrayList<>();
	    parentNodes.addAll(this.subGraphNodesMap.map.values());
	    return parentNodes;
	}

	/**
	 * Retrieves a list of parent nodes corresponding to the provided list of child nodes.
	 *
	 * @param childNodes The list of child nodes for which parent nodes are to be retrieved.
	 * @return A list of parent nodes corresponding to the provided child nodes.
	 */
	public ArrayList<NodeGraph> getParentNodes(ArrayList<NodeGraph> childNodes) {
	    final ArrayList<NodeGraph> parentNodes = new ArrayList<>();
	    for (final NodeGraph child : childNodes) {
	        final NodeGraph parent = this.subGraphNodesMap.findParent(child);
	        if (parent != null)
	            parentNodes.add(parent);
	    }
	    return parentNodes;
	}

	/**
	 * Retrieves a list of child nodes corresponding to the provided list of parent nodes.
	 *
	 * @param parentNodes The list of parent nodes for which child nodes are to be retrieved.
	 * @return A list of child nodes corresponding to the provided parent nodes.
	 */
	public ArrayList<NodeGraph> getChildNodes(ArrayList<NodeGraph> parentNodes) {
	    final ArrayList<NodeGraph> childNodes = new ArrayList<>();
	    for (final NodeGraph parent : parentNodes) {
	        final NodeGraph child = this.subGraphNodesMap.findChild(parent);
	        if (child != null)
	            childNodes.add(child);
	    }
	    return childNodes;
	}

	/**
	 * Retrieves the child edges corresponding to the provided parent edges.
	 *
	 * @param parentEdges The list of parent edges for which child edges are to be retrieved.
	 * @return A list of child edges corresponding to the provided parent edges.
	 */
	public ArrayList<EdgeGraph> getChildEdges(ArrayList<EdgeGraph> parentEdges) {
	    final ArrayList<EdgeGraph> childEdges = new ArrayList<>();
	    for (final EdgeGraph parent : parentEdges) {
	        final EdgeGraph child = this.subGraphEdgesMap.findChild(parent);
	        if (child != null)
	            childEdges.add(child);
	    }
	    return childEdges;
	}

	/**
	 * Retrieves the parent edge corresponding to the provided child edge.
	 *
	 * @param childEdge The child edge for which the parent edge is to be retrieved.
	 * @return The parent edge corresponding to the provided child edge, or null if not found.
	 */
	public EdgeGraph getParentEdge(EdgeGraph childEdge) {
	    return this.subGraphEdgesMap.findParent(childEdge);
	}

	/**
	 * Retrieves the parent edges corresponding to the provided child edges.
	 *
	 * @param childEdges The list of child edges for which parent edges are to be retrieved.
	 * @return A list of parent edges corresponding to the provided child edges.
	 */
	public ArrayList<EdgeGraph> getParentEdges(ArrayList<EdgeGraph> childEdges) {
	    final ArrayList<EdgeGraph> parentEdges = new ArrayList<>();
	    for (final EdgeGraph child : childEdges) {
	        final EdgeGraph parent = this.subGraphEdgesMap.findParent(child);
	        if (parent != null)
	            parentEdges.add(parent);
	    }
	    return parentEdges;
	}

	/**
	 * Retrieves a list of nodes within the current subgraph.
	 *
	 * @return A list of nodes present in the current subgraph.
	 */
	public ArrayList<NodeGraph> getNodesList() {
	    final ArrayList<NodeGraph> nodesList = new ArrayList<>();
	    nodesList.addAll(this.subGraphNodesMap.map.keySet());
	    return nodesList;
	}

	/**
	 * Sets the barrier attributes for edges within the current subgraph.
	 * This method copies barrier attributes from the corresponding parent graph's edges
	 * to the edges within the subgraph, and it compiles a list of unique barrier IDs
	 * present in the subgraph.
	 *
	 * @return A list of unique barrier IDs in the current subgraph.
	 */
	public void setAttributesChildEdge(EdgeGraph childEdge, EdgeGraph parentEdge) {

		childEdge.setID(parentEdge.getID());
		childEdge.dualNode = parentEdge.getDual();
		// for dual edges:
		childEdge.deflectionDegrees = parentEdge.deflectionDegrees;
	}

	/**
	 * It stores information about the barriers within this SubGraph.
	 *
	 */
	public void setSubGraphBarriers() {

		final ArrayList<Integer> graphBarriers = new ArrayList<>();
		for (final EdgeGraph childEdge : this.getEdges()) {
			childEdge.barriers = this.getParentEdge(childEdge).barriers;
			childEdge.positiveBarriers = this.getParentEdge(childEdge).positiveBarriers;
			childEdge.negativeBarriers = this.getParentEdge(childEdge).negativeBarriers;
			childEdge.waterBodies = this.getParentEdge(childEdge).waterBodies;
			childEdge.parks = this.getParentEdge(childEdge).parks;
			graphBarriers.addAll(childEdge.barriers);
		}
		final Set<Integer> setBarriers = new HashSet<>(graphBarriers);
		this.graphBarriers = new ArrayList<>(setBarriers);
	}

	/**
	 * Sets landmarks and visibility attributes for nodes within the current subgraph.
	 * This method copies landmarks and visibility attributes from the corresponding
	 * parent graph's nodes to the nodes within the subgraph.
	 */
	public void setSubGraphLandmarks() {
		final ArrayList<NodeGraph> childNodes = this.getNodesList();

		for (final NodeGraph node : childNodes) {
			final NodeGraph parentNode = this.getParentNode(node);
			node.visible2d = parentNode.visible2d;
			node.localLandmarks = parentNode.localLandmarks;
			node.distantLandmarks = parentNode.distantLandmarks;
			node.anchors = parentNode.anchors;
			node.distances = parentNode.distances;
		}
	}

	/**
	 * Returns the list of barrier IDs associated with the current subgraph. These barriers
	 * represent physical obstacles or features within the subgraph that shape movement and agent's cognitive maps.
	 *
	 * @return The list of barrier IDs within the subgraph.
	 */
	public ArrayList<Integer> getSubGraphBarriers() {
		return this.graphBarriers;
	}

	/**
	 * Generates the centrality map for the current subgraph by calculating the centrality values
	 * of its nodes based on the centrality values of their corresponding parent nodes in the parent
	 * graph. The resulting centrality values for nodes in the subgraph are sorted in descending order
	 * and stored as the subgraph's centrality map.
	 */
	public void generateSubGraphCentralityMap() {
		final LinkedHashMap<NodeGraph, Double> centralityMap = new LinkedHashMap<>();
		final Collection<NodeGraph> nodes = this.subGraphNodesMap.map.keySet();
		for (final NodeGraph n : nodes) {
			final NodeGraph parentNode = this.getParentNode(n);
			centralityMap.put(n, parentNode.centrality);
		}
		centralityMap = (LinkedHashMap<NodeGraph, Double>) Utilities.sortByValue(centralityMap, false);
	}

	/**
	 * Identifies and returns the global salient nodes within the current subgraph based on a specified
	 * percentile of centrality values in the parent graph. This method calculates salient nodes within
	 * the parent graph, filters them to retain only those that are parent nodes of the current subgraph,
	 * and returns the resulting list of global salient nodes in the subgraph.
	 *
	 * @param percentile The percentile value used to determine salient nodes in the parent graph.
	 * @return A list of global salient nodes within the current subgraph based on the specified percentile
	 *         in the parent graph.
	 */
	public ArrayList<NodeGraph> globalSalientNodesInSubGraph(double percentile) {
		final Map<NodeGraph, Double> parentGraphSalientNodesMap = parentGraph.graphSalientNodes(percentile);
		final ArrayList<NodeGraph> parentGraphSalientNodes = new ArrayList<>(parentGraphSalientNodesMap.keySet());
		parentGraphSalientNodes.retainAll(this.getParentNodes());
		return parentGraphSalientNodes;
	}
	
	/**
	 * Sets the salient nodes within the current subgraph based on a specified percentile
	 * of centrality values. This method calculates and stores the salient nodes within
	 * the subgraph by invoking the subGraphSalientNodes method with the provided percentile.
	 *
	 * @param salientNodesPercentile The percentile value used to determine salient nodes
	 *                               within the subgraph.
	 */
	public void setSubGraphSalientNodes(double salientNodesPercentile) {
		salientNodes = subGraphSalientNodes(salientNodesPercentile);
	}

	/**
	 * Calculates and returns salient nodes within a subgraph based on a specified percentile
	 * of centrality values. This method identifies nodes with centrality values equal to or
	 * higher than the specified percentile within the subgraph, then maps those nodes to
	 * their parent nodes in the parent graph and returns the resulting subgraph's salient nodes.
	 *
	 * @param percentile The percentile value used to determine salient nodes.
	 * @return A mapping of salient nodes within the subgraph to their centrality values, or
	 *         null if no salient nodes are found in the subgraph.
	 */
	private Map<NodeGraph, Double> subGraphSalientNodes(double percentile) {

		int position;
		position = (int) (centralityMap.size() * percentile);
		final double boundary = new ArrayList<>(centralityMap.values()).get(position);

		final Map<NodeGraph, Double> filteredMap = centralityMap.entrySet().stream()
				.filter(entry -> entry.getValue() >= boundary)
				.collect(Collectors.toMap(entry -> entry.getKey(), entry -> entry.getValue()));

		if (filteredMap.size() == 0 || filteredMap == null)
			return null;
		final Map<NodeGraph, Double> parentMap = new HashMap<>();

		for (final Map.Entry<NodeGraph, Double> entry : filteredMap.entrySet()) {
			final NodeGraph parentNode = this.getParentNode(entry.getKey());
			parentMap.put(parentNode, entry.getValue());
		}
		return parentMap;
	}
}