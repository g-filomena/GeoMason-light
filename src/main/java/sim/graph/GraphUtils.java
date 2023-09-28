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

import org.locationtech.jts.geom.Coordinate;
import org.locationtech.jts.geom.Geometry;
import org.locationtech.jts.geom.GeometryFactory;
import org.locationtech.jts.geom.LineString;
import org.locationtech.jts.geom.Point;

/**
 * The {@code GraphUtils} class provides utility methods for working {@code NodeGraph} objects.
 */
public class GraphUtils {

	/**
	 * Calculates the Euclidean distance between two nodes in a graph.
	 *
	 * @param node       The first node.
	 * @param otherNode  The second node.
	 * @return The Euclidean distance between the two nodes' coordinates.
	 */
	public static double nodesDistance(NodeGraph node, NodeGraph otherNode) {
		final Coordinate originCoord = node.getCoordinate();
		final Coordinate destinationCoord = otherNode.getCoordinate();
		return euclideanDistance(originCoord, destinationCoord);
	}

	/**
	 * Calculates the minimum enclosing circle (smallest circle that completely encloses a collection of nodes).
	 *
	 * @param nodes The collection of nodes for which to calculate the enclosing circle.
	 * @return A geometry representing the minimum enclosing circle.
	 */
	public static Geometry enclosingCircleFromNodes(Iterable<NodeGraph> nodes) {
		// Calculate the center of the enclosing circle
		Coordinate center = calculateCenterPointNodes(nodes);
		// Calculate the radius of the enclosing circle
		double radius = calculateRadius(center, nodes);

		// Create a circle geometry representing the enclosing circle
		GeometryFactory geometryFactory = new GeometryFactory();
		Point circleCenter = geometryFactory.createPoint(center);
		Geometry enclosingCircle = circleCenter.buffer(radius);
		return enclosingCircle;
	}
	
	/**
	 * Calculates the smallest enclosing circle for two nodes in a graph.
	 *
	 * @param node       The first node.
	 * @param otherNode  The second node.
	 * @return The smallest enclosing circle as a geometry.
	 */
	public static Geometry enclosingCircleBetweenNodes(NodeGraph node, NodeGraph otherNode) {
		final LineString line = LineStringBetweenNodes(node, otherNode);
		final Point centroid = line.getCentroid();
		final Geometry smallestEnclosingCircle = centroid.buffer(line.getLength() / 2);
		return smallestEnclosingCircle;
	}

    /**
     * Calculates the center point (centroid) of a collection of nodes.
     *
     * @param nodes The collection of nodes for which to calculate the center point.
     * @return The center point coordinate.
     */
    private static Coordinate calculateCenterPointNodes(Iterable<NodeGraph> nodes) {
        double totalX = 0.0;
        double totalY = 0.0;
        int count = 0;

        for (NodeGraph node : nodes) {
            totalX += node.getCoordinate().getX();
            totalY += node.getCoordinate().getY();
            count++;
        }

        return new Coordinate(totalX / count, totalY / count);
    }

    /**
     * Calculates the radius of a collection of nodes relative to a center point.
     *
     * @param center The center point coordinate.
     * @param nodes  The collection of nodes for which to calculate the radius.
     * @return The maximum distance from the center point to any node.
     */
    private static double calculateRadius(Coordinate center, Iterable<NodeGraph> nodes) {
        double maxDistance = 0.0;

        for (NodeGraph node : nodes) {
            double distance = center.distance(node.getCoordinate());
            if (distance > maxDistance) {
                maxDistance = distance;
            }
        }
        return maxDistance;
    }

    /**
     * Finds the closest node to a target coordinate from a collection of nodes.
     *
     * @param targetCoordinates The target coordinate to which to find the closest node.
     * @param nodes             The collection of nodes to search for the closest node.
     * @return The closest node to the target coordinate.
     */
    public static NodeGraph findClosestNode(Coordinate targetCoordinates, Iterable<NodeGraph> nodes) {
        NodeGraph closestNode = null;
        double minDistance = Double.MAX_VALUE;

        for (NodeGraph node : nodes) {
            double distance = targetCoordinates.distance(node.getCoordinate());
            if (distance < minDistance) {
                minDistance = distance;
                closestNode = node;
            }
        }
        return closestNode;
    }
    
    /**
	 * It returns a LineString between two given nodes.
	 *
	 * @param node      a node;
	 * @param otherNode an other node;
	 */
	public static LineString LineStringBetweenNodes(NodeGraph node, NodeGraph otherNode) {
		final Coordinate[] coords = {node.getCoordinate(), otherNode.getCoordinate()};
		final LineString line = new GeometryFactory().createLineString(coords);
		return line;
	}
	
	/**
	 * Computes the Euclidean distance between two locations
	 *
	 * @param originCoord      the origin location;
	 * @param destinationCoord the destination;
	 */
	public static double euclideanDistance(Coordinate originCoord, Coordinate destinationCoord) {
		return Math.sqrt(
				Math.pow(originCoord.x - destinationCoord.x, 2) + Math.pow(originCoord.y - destinationCoord.y, 2));
	}
}
