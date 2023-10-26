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

import java.util.HashMap;
import java.util.Map;

import org.locationtech.jts.geom.Coordinate;
import org.locationtech.jts.geom.LineString;
import org.locationtech.jts.planargraph.Edge;

import sim.util.geo.AttributeValue;
import sim.util.geo.MasonGeometry;

/**
 * An edge that extends the and `Edge` (Jts) class. This is one of the two
 * components, along with {@link NodeGraph}, of graphs belonging to the class
 * {@link Graph}.
 */
public class EdgeGraph extends Edge {

	public int edgeID;
	public int regionID;
	public double deflectionDegrees;
	public MasonGeometry masonGeometry;
	private LineString line; // line that corresponds to this edge

	public NodeGraph fromNode, toNode;
	public NodeGraph dualNode;

	public HashMap<String, Integer> volumes = new HashMap<>();
	public Map<String, AttributeValue> attributes;
	private Coordinate centroidCoords;
	private Double length;
	public boolean isKnown = false;

	/**
	 * Sets the attributes for this EdgeGraph.
	 *
	 * @param attributes A map of attribute names and values.
	 */
	public void setAttributes(final Map<String, AttributeValue> attributes) {
		this.attributes = attributes;
	}

	/**
	 * Constructs an EdgeGraph with the provided LineString.
	 *
	 * @param line The LineString representing the edge.
	 */
	public EdgeGraph(LineString line) {
		this.line = line;
		length = line.getLength();
		centroidCoords = line.getCentroid().getCoordinate();
	}

	/**
	 * Gets the LineString representing this edge.
	 *
	 * @return The LineString for this edge.
	 */
	public LineString getLine() {
		return line;
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
	public Integer getID() {
		return this.edgeID;
	}

	/**
	 * Returns the edge's corresponding dual node.
	 *
	 * @return The dual node associated with this edge.
	 */
	public NodeGraph getDual() {
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
	 * Returns the deflection angle if this edge is a dual edge and represents a
	 * link between two dual nodes (street segment).
	 *
	 * @return The deflection angle of the dual edge, if applicable; otherwise, it
	 *         returns 0.0.
	 */
	public double getDeflectionAngle() {
		return deflectionDegrees;
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
	 * Resets the volumes when called during the simulation (e.g., at the start of a
	 * new run).
	 */
	public void resetVolumes() {

		for (final String key : volumes.keySet())
			volumes.replace(key, 0);
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
