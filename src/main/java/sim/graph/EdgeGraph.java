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
import org.locationtech.jts.geom.LineString;
import org.locationtech.jts.planargraph.DirectedEdge;
import org.locationtech.jts.planargraph.Edge;

import sim.util.geo.AttributeValue;
import sim.util.geo.MasonGeometry;

/**
 * Represents an edge in a graph, extending the JTS (Java Topology Suite) Edge
 * class. This class is a component of graphs that are part of the {@link Graph}
 * class, along with {@link NodeGraph}. It encapsulates additional attributes
 * and functionalities specific to the graph's edges.
 */
public class EdgeGraph extends Edge {

	protected int edgeID;
	protected int regionID = -1;
	protected MasonGeometry masonGeometry;
	private final LineString line; // line that corresponds to this edge

	private NodeGraph fromNode;
	private NodeGraph toNode;
	private List<NodeGraph> nodes = new ArrayList<NodeGraph>();
	protected NodeGraph dualNode;
	private final Coordinate centroidCoords;
	private final Double length;
	public Map<String, AttributeValue> attributes = new HashMap<>();
	public boolean isKnown = false;

	// Attributes specific to dual edges
	protected double deflectionAngle;

	/**
	 * Sets the attributes for this EdgeGraph.
	 *
	 * @param attributes A map of attribute names and their corresponding values.
	 */
	public void setAttributes(final Map<String, AttributeValue> attributes) {
		this.attributes = attributes;
	}

	/**
	 * Constructs an EdgeGraph with the specified LineString.
	 *
	 * @param line The LineString representing the physical geometry of the edge.
	 */
	public EdgeGraph(LineString line) {
		this.line = line;
		length = line.getLength();
		centroidCoords = line.getCentroid().getCoordinate();
	}

	/**
	 * Initializes this Edge's two DirectedEdges, and for each DirectedEdge: sets
	 * the Edge, sets the symmetric DirectedEdge, and adds this Edge to its
	 * from-Node.
	 */
	@Override
	public void setDirectedEdges(DirectedEdge de0, DirectedEdge de1) {
		dirEdge = new DirectedEdge[] { de0, de1 };
		de0.setEdge(this);
		de1.setEdge(this);
		de0.setSym(de1);
		de1.setSym(de0);
		de0.getFromNode().addOutEdge(de0);
		de1.getFromNode().addOutEdge(de1);
	}

	/**
	 * Gets the LineString representing this edge.
	 *
	 * @return The LineString for this edge.
	 */
	public LineString getLine() {
		return (LineString) this.line.copy();
	}

	public Coordinate[] getCoordinates() {
		return this.line.getCoordinates();
	}

	/**
	 * Sets the MasonGeometry of the edge.
	 *
	 */
	public void setMasonGeometry(MasonGeometry masonGeometry) {
		this.masonGeometry = masonGeometry;
	}

	/**
	 * Returns the MasonGeometry of the edge.
	 *
	 * @return The MasonGeometry of the edge.
	 */
	public MasonGeometry getMasonGeometry() {
		return masonGeometry;
	}

	/**
	 * Sets the nodes connected by this edge.
	 *
	 * @param fromNode The starting node of the edge.
	 * @param toNode   The ending node of the edge.
	 */
	public void setNodes(final NodeGraph fromNode, final NodeGraph toNode) {
		this.fromNode = fromNode;
		this.toNode = toNode;
		nodes.add(fromNode);
		nodes.add(toNode);
	}

	/**
	 * Retrieves the list of nodes in the graph.
	 *
	 * @return a list of nodes in the graph
	 */
	public List<NodeGraph> getNodes() {
		return nodes;
	}

	/**
	 * Retrieves the starting node of the edge.
	 *
	 * @return the starting node of the edge
	 */
	public NodeGraph getFromNode() {
		return fromNode;
	}

	/**
	 * Retrieves the ending node of the edge.
	 *
	 * @return the ending node of the edge
	 */
	public NodeGraph getToNode() {
		return toNode;
	}

	/**
	 * Sets the ID of the edge.
	 *
	 * @param edgeID The unique identifier to assign to the edge.
	 */
	public void setID(int edgeID) {
		this.edgeID = edgeID;
	}

	/**
	 * Gets the ID of this edge.
	 *
	 * @return The ID of this edge.
	 */
	public int getID() {
		return this.edgeID;
	}

	/**
	 * Sets the ID of the edge's region.
	 *
	 * @param regionID The ID to set for the edge's region.
	 */
	public void setRegionID(int regionID) {
		this.regionID = regionID;
	}

	/**
	 * Gets the ID of the region of this edge.
	 *
	 * @return The ID of the region of this edge.
	 */
	public int getRegionID() {
		return this.regionID;
	}

	/**
	 * Sets the dual node to the specified centroid.
	 *
	 * @param centroid the node to be set as the dual node
	 */
	public void setDualNode(NodeGraph centroid) {
		dualNode = centroid;
	}

	/**
	 * Returns the edge's corresponding dual node.
	 *
	 * @return The dual node associated with this edge.
	 */
	public NodeGraph getDualNode() {
		return this.dualNode;
	}

	/**
	 * Returns the edge's length.
	 *
	 * @return The length of the edge.
	 */
	public double getLength() {
		return length;
	}

	/**
	 * Set the deflection angle; this would make sense if this edge is a dual edge
	 * and represents a link between two dual nodes (street segments).
	 *
	 * @param angle the value to be used to set the angle.
	 */
	public void setDeflectionAngle(Double angle) {
		deflectionAngle = angle;
	}

	/**
	 * Returns the deflection angle if this edge is a dual edge and represents a
	 * link between two dual nodes (street segment).
	 *
	 * @return The deflection angle of the dual edge, if applicable; otherwise, it
	 *         returns 0.0.
	 */
	public double getDeflectionAngle() {
		return deflectionAngle;
	}

	/**
	 * Returns the coordinate of the edge's centroid.
	 *
	 * @return The coordinate of the edge's centroid.
	 */
	public Coordinate getCoordsCentroid() {
		return centroidCoords;
	}

	/**
	 * Given one of the nodes of this edge, it returns the other one.
	 *
	 * @param node One of the nodes connected to this edge.
	 * @return The other node connected to this edge, or null if the provided node
	 *         is not connected to this edge.
	 */
	public NodeGraph getOtherNode(NodeGraph node) {
		if (fromNode.equals(node))
			return toNode;
		else if (toNode.equals(node))
			return fromNode;
		else
			return null;
	}

	/**
	 * Gets the common node between this edge and another edge.
	 *
	 * @param edge Another edge to check for a common node.
	 * @return The common node if found, or null if there is no common node between
	 *         the edges.
	 */
	public NodeGraph getCommonNode(EdgeGraph edge) {

		if (fromNode.equals(edge.toNode) || fromNode.equals(edge.fromNode))
			return fromNode;
		if (toNode.equals(edge.toNode) || toNode.equals(edge.fromNode))
			return toNode;
		else
			return null;
	}
}
