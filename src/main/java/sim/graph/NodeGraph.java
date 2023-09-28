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
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.locationtech.jts.geom.Coordinate;
import org.locationtech.jts.planargraph.DirectedEdge;
import org.locationtech.jts.planargraph.Node;

import sim.util.geo.MasonGeometry;
import sim.util.geo.Utilities;

/**
 * A node of a planar graph that extends the {@code Node} class (JTS), with additional
 * and more straightforward functions. This class is one of the two components,
 * along with {@code EdgeGraph}, of the graphs belonging to the {@code Graph} class.
 */
public class NodeGraph extends Node {

    /**
     * Constructs a new NodeGraph instance with the specified coordinate.
     *
     * @param pt The coordinate associated with this node.
     */
    public NodeGraph(Coordinate pt) {
        super(pt);
    }
	public int nodeID;
	public int regionID = -1;
	public boolean gateway;

	public MasonGeometry masonGeometry;
	public EdgeGraph primalEdge;
	public double centrality, centrality_sc;

	public ArrayList<Building> visible2d = new ArrayList<>();
	public ArrayList<Building> localLandmarks = new ArrayList<>();
	public ArrayList<Building> distantLandmarks = new ArrayList<>();
	public ArrayList<Building> anchors = new ArrayList<>();
	public List<Double> distances = new ArrayList<>();

	public List<Integer> adjacentRegions = new ArrayList<>();
	public ArrayList<NodeGraph> adjacentEntries = new ArrayList<>();
	private ArrayList<NodeGraph> adjacentNodes = new ArrayList<>();
	private ArrayList<EdgeGraph> edges = new ArrayList<>();
	private ArrayList<DirectedEdge> outEdges = new ArrayList<>();
	public String DMA = "";

	/**
	 * Sets the ID of the node.
	 *
	 * @param nodeID The ID to set for the node.
	 */
	public void setID(int nodeID) {
	    this.nodeID = nodeID;
	}

	/**
	 * Returns the ID of the node.
	 *
	 * @return The ID of the node.
	 */
	public Integer getID() {
	    return nodeID;
	}

	/**
	 * Returns the list of edges that depart from this node.
	 *
	 * @return The list of edges departing from the node.
	 */
	public ArrayList<EdgeGraph> getEdges() {
	    return edges;
	}

	/**
	 * Returns the list of nodes that are reachable from this node.
	 *
	 * @return The list of nodes that can be reached from this node.
	 */
	public ArrayList<NodeGraph> getAdjacentNodes() {
	    return adjacentNodes;
	}

	/**
	 * Sets lists useful for identifying the neighboring components of this node.
	 * This method initializes the edge and adjacent node lists.
	 */
	public void setNeighbouringComponents() {
	    setEdgesNode();
	    setAdjacentNodes();
	}

	/**
	 * Returns the list of directed edges that depart from the node (out-going edges).
	 *
	 * @return The list of out-going directed edges from the node.
	 */
	public ArrayList<DirectedEdge> getOutsDirectedEdges() {
	    return outEdges;
	}

	/**
	 * Identifies and sets the list of all the edges (non-directed) that depart from this node.
	 * This method initializes the 'edges' list with non-duplicated edges.
	 */
	private void setEdgesNode() {
	    final ArrayList<EdgeGraph> edges = new ArrayList<>();
	    this.outEdges = new ArrayList<DirectedEdge>(this.getOutEdges().getEdges());

	    for (final DirectedEdge outEdge : outEdges) {
	        final EdgeGraph edge = (EdgeGraph) outEdge.getEdge();
	        if (!edges.contains(edge))
	            edges.add(edge);
	    }
	    this.edges = edges;
	}

	/**
	 * Identifies a list of all the nodes adjacent to this node (i.e., sharing an edge with this node).
	 * This method initializes the 'adjacentNodes' list.
	 */
	private void setAdjacentNodes() {
	    final ArrayList<NodeGraph> adjacentNodes = new ArrayList<>();

	    for (final EdgeGraph edge : edges) {
	        final NodeGraph opposite = (NodeGraph) edge.getOppositeNode(this);
	        adjacentNodes.add(opposite);
	    }
	    this.adjacentNodes = adjacentNodes;
	}

	/**
	 * Returns a list of integers representing all the adjacent regions to this node, if any.
	 * If this node is not a gateway, it returns null.
	 *
	 * @return A list of integers representing adjacent regions, or null if this node is not a gateway.
	 */
	public ArrayList<Integer> getAdjacentRegion() {
	    if (!gateway)
	        return null;

	    final ArrayList<NodeGraph> oppositeNodes = new ArrayList<>(adjacentNodes);
	    final ArrayList<Integer> adjacentRegions = new ArrayList<>();

	    for (final NodeGraph opposite : oppositeNodes) {
	        final int regionID = opposite.regionID;
	        if (regionID != this.regionID) {
	            adjacentRegions.add(regionID);
	            adjacentEntries.add(opposite);
	        }
	    }
	    return adjacentRegions;
	}

	/**
	 * Returns a node in the dual representation of the graph, which represents a
	 * segment departing from this node (a node in a dual graph represents an actual
	 * street segment). When the node for which the dual-node is desired is the
	 * originNode of a desired route, the segment closest to the destination is
	 * chosen (~ towards it). When the node for which the dual-node is desired is
	 * the destinationNode, the segment closest to the origin is chosen.
	 *
	 * @param originNode            The origin node of the desired route.
	 * @param destinationNode       The destination node of the desired route.
	 * @param regionBasedNavigation Indicates if the agent is navigating through regions.
	 * @param previousJunction      Node to avoid choosing dual nodes representing segments that should not be traversed.
	 * @return The selected dual node representing a street segment.
	 */
	public NodeGraph getDualNode(NodeGraph originNode, NodeGraph destinationNode, boolean regionBasedNavigation,
			NodeGraph previousJunction) {
		double distance = Double.MAX_VALUE;
		NodeGraph best = null;
		NodeGraph dualNode = null;
		final ArrayList<EdgeGraph> edges = new ArrayList<>(this.edges);
		double cost;

		for (final EdgeGraph edge : edges) {
			if (edge.regionID == -1 && regionBasedNavigation)
				continue;
			dualNode = edge.getDual();
			if (dualNode == null)
				continue;

			if (this.equals(destinationNode)) {
				cost = NodeGraphUtils.nodesDistance(edge.getOtherNode(this), originNode);
				if (cost < distance) {
					distance = cost;
					best = dualNode;
				}
			} else {
				cost = NodeGraphUtils.nodesDistance(edge.getOtherNode(this), destinationNode);

				if (previousJunction != null && (previousJunction == dualNode.primalEdge.fromNode
						|| previousJunction == dualNode.primalEdge.toNode))
					continue;
				if (regionBasedNavigation) {
					final ArrayList<EdgeGraph> nextEdges = new ArrayList<>(edge.getOtherNode(originNode).getEdges());
					nextEdges.remove(edge);
					boolean bridges = true;
					for (final EdgeGraph next : nextEdges)
						if (next.regionID != -1)
							bridges = false;
					if (bridges)
						continue;
				}
				if (cost < distance) {
					distance = cost;
					best = dualNode;
				}
			}
		}
		return best;
	}

	/**
	 * Returns a map of nodes in the dual representation of the graph, which
	 * represent segments departing from this node (a node in a dual graph
	 * represents an actual street segment). When the node for which the dual-nodes
	 * are desired is the originNode of a desired route, the map is sorted on the basis
	 * of the distance from the dual nodes (segments' centroids) to the destination
	 * node (~ towards it). When the node for which the dual-nodes are desired is
	 * the destinationNode, the map is sorted on the basis of the distance from the
	 * dual nodes (segments' centroids) to the origin node.
	 *
	 * This function is preferable to the one above as it allows considering
	 * multiple segments as possible "departures" for route planning algorithms.
	 * When computing paths within subgraphs, some specific segments may be indeed
	 * unreachable.
	 *
	 * @param originNode            The origin node of the desired route.
	 * @param destinationNode       The destination node of the desired route.
	 * @param regionBasedNavigation Indicates if the agent is navigating through regions.
	 * @param previousJunction      Node to avoid choosing dual nodes representing segments that should not be traversed.
	 * @return A map of selected dual nodes representing street segments, sorted by distance.
	 */
	public Map<NodeGraph, Double> getDualNodes(NodeGraph originNode, NodeGraph destinationNode,
			boolean regionBasedNavigation, NodeGraph previousJunction) {
		final HashMap<NodeGraph, Double> dualNodes = new HashMap<>();
		NodeGraph dualNode = null;
		final ArrayList<EdgeGraph> edges = new ArrayList<>(this.edges);

		for (final EdgeGraph edge : edges) {
			if (edge.regionID == -1 && regionBasedNavigation)
				continue;
			dualNode = edge.getDual();
			if (dualNode == null)
				continue;

			if (this.equals(destinationNode)) {
				final double cost = NodeGraphUtils.nodesDistance(edge.getOtherNode(this), originNode);
				dualNodes.put(dualNode, cost);
			} else {
				final double cost = NodeGraphUtils.nodesDistance(edge.getOtherNode(this), destinationNode);

				if (previousJunction != null && (previousJunction == dualNode.primalEdge.fromNode
						|| previousJunction == dualNode.primalEdge.toNode))
					continue;
				if (regionBasedNavigation) {
					final ArrayList<EdgeGraph> nextEdges = new ArrayList<>(edge.getOtherNode(originNode).getEdges());
					nextEdges.remove(edge);
					boolean bridges = true;
					for (final EdgeGraph next : nextEdges)
						if (next.regionID != -1)
							bridges = false;
					if (bridges)
						continue;
				}
				dualNodes.put(dualNode, cost);
			}
		}

		return Utilities.sortByValue(dualNodes, false);
	}
}
