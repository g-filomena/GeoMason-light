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
import java.util.Arrays;
import java.util.Comparator;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.concurrent.ConcurrentHashMap;
import java.util.stream.Collectors;
import java.util.stream.StreamSupport;

import org.javatuples.Pair;
import org.locationtech.jts.geom.Coordinate;
import org.locationtech.jts.geom.Geometry;
import org.locationtech.jts.geom.GeometryFactory;
import org.locationtech.jts.geom.LineString;
import org.locationtech.jts.geom.LinearRing;
import org.locationtech.jts.geom.MultiLineString;
import org.locationtech.jts.geom.Point;
import org.locationtech.jts.geom.Polygon;
import org.locationtech.jts.operation.union.UnaryUnionOp;
import org.locationtech.jts.planargraph.DirectedEdge;

import sim.field.geo.VectorLayer;
import sim.util.geo.GeometryUtilities;
import sim.util.geo.MasonGeometry;

/**
 * The `GraphUtils` class provides utility methods for working {@link NodeGraph}
 * objects.
 */
public class GraphUtils {

	// Introduce a HashMap to cache previously computed distances
	private static Map<Pair<Coordinate, Coordinate>, Double> distanceCache = new ConcurrentHashMap<>();
	private static Map<Pair<Coordinate, Coordinate>, Polygon> visibilityPolygonsCache = new ConcurrentHashMap<>();
	private static final GeometryFactory GEOMETRY_FACTORY = new GeometryFactory();
	static int DISTANCE_ALONG_VISIBILITY = 10;

	/**
	 * Calculates the Euclidean distance between two nodes in a graph.
	 *
	 * @param node      The first node.
	 * @param otherNode The second node.
	 * @return The Euclidean distance between the two nodes' coordinates.
	 */
	public static double nodesDistance(NodeGraph node, NodeGraph otherNode) {

		final Coordinate coords = node.getCoordinate();
		final Coordinate otherCoords = otherNode.getCoordinate();
		Pair<Coordinate, Coordinate> pair = new Pair<>(coords, otherCoords);
		Pair<Coordinate, Coordinate> otherPair = new Pair<>(otherCoords, coords);

		if (distanceCache.containsKey(pair))
			return distanceCache.get(pair);
		else if (distanceCache.containsKey(otherPair))
			return distanceCache.get(otherPair);

		double distance = GeometryUtilities.euclideanDistance(coords, otherCoords);
		distanceCache.put(pair, distance);
		return distance;
	}

	/**
	 * Calculates the minimum enclosing circle (smallest circle that completely
	 * encloses a collection of nodes).
	 *
	 * @param nodes The collection of nodes for which to calculate the enclosing
	 *              circle.
	 * @return A geometry representing the minimum enclosing circle.
	 */
	public static Geometry smallestEnclosingGeometryBetweenNodes(List<NodeGraph> nodes) {

		if (nodes.size() == 1)
			return nodes.get(0).masonGeometry.getGeometry().buffer(50);
		if (nodes.size() == 2)
			return enclosingCircleBetweenTwoNodes(nodes.get(0), nodes.get(1));
		else
			return convexHullFromNodes(nodes);
	}

	/**
	 * Calculates the convex hull from a list of nodes using Andrew's monotone chain
	 * algorithm.
	 *
	 * @param nodes The list of nodes to compute the convex hull from.
	 * @return The convex hull polygon as a Geometry object, or null if there are
	 *         fewer than 3 nodes.
	 */
	private static Geometry convexHullFromNodes(List<NodeGraph> nodes) {

		// Sort nodes by x-coordinate (break ties by y-coordinate)
		nodes.sort(Comparator.comparingDouble((NodeGraph node) -> node.getCoordinate().x)
				.thenComparingDouble(node -> node.getCoordinate().y));

		List<NodeGraph> hull = new ArrayList<>();

		// Build lower hull
		for (NodeGraph node : nodes) {
			while (hull.size() >= 2 && !isCounterClockwise(hull.get(hull.size() - 2), hull.get(hull.size() - 1), node))
				hull.remove(hull.size() - 1);
			hull.add(node);
		}

		// Build upper hull
		int lowerHullSize = hull.size();
		for (int i = nodes.size() - 1; i >= 0; i--) {
			NodeGraph node = nodes.get(i);
			while (hull.size() > lowerHullSize
					&& !isCounterClockwise(hull.get(hull.size() - 2), hull.get(hull.size() - 1), node))
				hull.remove(hull.size() - 1);

			hull.add(node);
		}

		// Remove the last point because it is repeated at the beginning of the list
		if (hull.size() > 1 && hull.get(hull.size() - 1).equals(hull.get(0)))
			hull.remove(hull.size() - 1);

		// Create an array of coordinates for the convex hull
		Coordinate[] coordinates = hull.stream().map(NodeGraph::getCoordinate).toArray(Coordinate[]::new);

		// Ensure the polygon is closed
		if (!coordinates[0].equals(coordinates[coordinates.length - 1])) {
			coordinates = Arrays.copyOf(coordinates, coordinates.length + 1);
			coordinates[coordinates.length - 1] = coordinates[0];
		}

		// Create and return the polygon using JTS
		LinearRing linearRing = GEOMETRY_FACTORY.createLinearRing(coordinates);
		return GEOMETRY_FACTORY.createPolygon(linearRing);
	}

	/**
	 * Determines if three nodes form a counter-clockwise turn.
	 *
	 * This method uses the cross product of vectors to determine the relative
	 * orientation of three points (nodes). It returns true if the points form a
	 * counter-clockwise turn, and false otherwise.
	 *
	 * @param a The first node.
	 * @param b The second node.
	 * @param c The third node.
	 * @return true if the nodes a, b, and c form a counter-clockwise turn, false
	 *         otherwise.
	 */
	private static boolean isCounterClockwise(NodeGraph a, NodeGraph b, NodeGraph c) {
		Coordinate p1 = a.getCoordinate();
		Coordinate p2 = b.getCoordinate();
		Coordinate p3 = c.getCoordinate();
		return (p2.x - p1.x) * (p3.y - p1.y) - (p2.y - p1.y) * (p3.x - p1.x) > 0;
	}

	/**
	 * Calculates the smallest enclosing circle for two nodes in a graph.
	 *
	 * @param node      The first node.
	 * @param otherNode The second node.
	 * @return The smallest enclosing circle as a geometry.
	 */
	protected static Geometry enclosingCircleBetweenTwoNodes(NodeGraph node, NodeGraph otherNode) {
		final LineString line = LineStringBetweenNodes(node, otherNode);
		final Point centroid = line.getCentroid();
		final Geometry smallestEnclosingCircle = centroid.buffer(line.getLength() / 2);
		return smallestEnclosingCircle;
	}

	/**
	 * Finds the closest node to a target coordinate from a collection of nodes.
	 *
	 * @param targetCoordinates The target coordinate to which to find the closest
	 *                          node.
	 * @param nodes             The collection of nodes to search for the closest
	 *                          node.
	 * @return The closest node to the target coordinate.
	 */
	public static NodeGraph findClosestNode(Coordinate targetCoordinates, Iterable<NodeGraph> nodes) {
		return StreamSupport.stream(nodes.spliterator(), true) // Convert Iterable to parallel stream
				.min(Comparator.comparingDouble(node -> targetCoordinates.distance(node.getCoordinate()))).orElse(null); // found
	}

	/**
	 * It returns a LineString between two given nodes.
	 *
	 * @param node      a node;
	 * @param otherNode an other node;
	 */
	public static LineString LineStringBetweenNodes(NodeGraph node, NodeGraph otherNode) {
		final Coordinate[] coords = { node.getCoordinate(), otherNode.getCoordinate() };
		final LineString line = new GeometryFactory().createLineString(coords);
		return line;
	}

	/**
	 * Given two centroids (nodes in the dual graph), identifies their shared
	 * junction (i.e., the junction shared by the corresponding primal segments).
	 *
	 * @param centroid      A dual node.
	 * @param otherCentroid Another dual node.
	 * @return The common primal junction node.
	 */
	public static NodeGraph getPrimalJunction(NodeGraph centroid, NodeGraph otherCentroid) {

		EdgeGraph edge = centroid.getPrimalEdge();
		EdgeGraph otherEdge = otherCentroid.getPrimalEdge();

		if (edge.getFromNode().equals(otherEdge.getFromNode()) || edge.getFromNode().equals(otherEdge.getToNode()))
			return edge.getFromNode();
		else if (edge.getToNode().equals(otherEdge.getFromNode()) || edge.getToNode().equals(otherEdge.getToNode()))
			return edge.getToNode();
		else
			return null;
	}

	/**
	 * Identifies the previous junction traversed in a dual graph path to avoid
	 * traversing an unnecessary segment in the primal graph.
	 *
	 * @param sequenceDirectedEdges A sequence of GeomPlanarGraphDirectedEdge
	 *                              representing the path.
	 * @return The previous junction node.
	 */
	public static NodeGraph previousJunction(List<DirectedEdge> sequenceDirectedEdges) {

		if (sequenceDirectedEdges.size() == 1)
			return (NodeGraph) sequenceDirectedEdges.get(0).getFromNode();

		int ixLast = sequenceDirectedEdges.size() - 1;
		int ixBeforeLast = sequenceDirectedEdges.size() - 2;
		NodeGraph lastCentroid = ((EdgeGraph) sequenceDirectedEdges.get(ixLast).getEdge()).getDualNode();
		NodeGraph otherCentroid = ((EdgeGraph) sequenceDirectedEdges.get(ixBeforeLast).getEdge()).getDualNode();
		return GraphUtils.getPrimalJunction(lastCentroid, otherCentroid);
	}

	public static Set<NodeGraph> nodesFromEdges(Set<EdgeGraph> edges) {
		return edges.stream().flatMap(edge -> edge.getNodes().stream()).collect(Collectors.toSet());
	}

	public static Set<EdgeGraph> edgesFromNodes(Set<NodeGraph> nodes) {
		return nodes.stream().flatMap(node -> node.getEdges().stream()).collect(Collectors.toSet());
	}

	public static List<Integer> getNodeIDs(List<NodeGraph> nodes) {
		return nodes.stream().map(NodeGraph::getID).collect(Collectors.toList());
	}

	public static List<Integer> getEdgeIDs(List<EdgeGraph> edges) {
		return edges.stream().map(EdgeGraph::getID).collect(Collectors.toList());
	}

	public static Polygon createVisibilityPolygon(NodeGraph fromNode, NodeGraph toNode, Double visibilityAngle,
			VectorLayer obstructions, double maxExpansionDistance) {

		Coordinate fromCoordinates = fromNode.getCoordinate();
		Coordinate toCoordinates = toNode.getCoordinate();

		Pair<Coordinate, Coordinate> pair = new Pair<>(fromCoordinates, toCoordinates);

		if (visibilityPolygonsCache.containsKey(pair))
			return visibilityPolygonsCache.get(pair);

		double edgeAngle = Math.toDegrees(Math.atan2(toCoordinates.getY() - fromCoordinates.getY(),
				toCoordinates.getX() - fromCoordinates.getX()));
		List<Coordinate> coords = new ArrayList<>();

		int visibilityLimit = (int) (visibilityAngle / 2.0);
		for (int i = -visibilityLimit; i <= visibilityLimit; i += DISTANCE_ALONG_VISIBILITY) {
			double angle = Math.toRadians(edgeAngle + i);
			double x = toCoordinates.getX() + maxExpansionDistance * Math.cos(angle);
			double y = toCoordinates.getY() + maxExpansionDistance * Math.sin(angle);
			coords.add(new Coordinate(x, y));
		}

		List<LineString> lines = new ArrayList<>();
		for (Coordinate coord : coords)
			lines.add(GEOMETRY_FACTORY.createLineString(new Coordinate[] { toNode.getCoordinate(), coord }));

		Geometry toNodeGeometry = toNode.getMasonGeometry().geometry;
		Geometry buffer = toNodeGeometry.buffer(maxExpansionDistance);
		List<MasonGeometry> possibleMatches = obstructions.containedFeatures(buffer);
		List<Geometry> obstacles = new ArrayList<>();

		possibleMatches.stream().map(masonGeometry -> masonGeometry.geometry)
				.filter(geometry -> lines.stream()
						.anyMatch(line -> geometry.crosses(line) && !geometry.touches(toNodeGeometry)))
				.forEach(obstacles::add);

		List<LineString> clippedLines = new ArrayList<>();
		if (!obstacles.isEmpty()) {
			Geometry unionObstacles = UnaryUnionOp.union(obstacles);

			for (LineString line : lines) {
				Geometry intersection = line.intersection(unionObstacles);
				if (!intersection.isEmpty()) {

					Coordinate[] intersectionCoords = intersection instanceof MultiLineString
							? intersection.getCoordinates()
							: new Coordinate[] { intersection.getCoordinate() };
					clippedLines.add(GEOMETRY_FACTORY
							.createLineString(new Coordinate[] { toNode.getCoordinate(), intersectionCoords[0] }));
				} else
					clippedLines.add(line);
			}
		} else
			clippedLines = lines;

		List<Coordinate> polyCoords = new ArrayList<>();
		polyCoords.add(toNode.getCoordinate());
		clippedLines.forEach(line -> polyCoords.add(line.getCoordinateN(1)));
		polyCoords.add(toNode.getCoordinate());

		Polygon visibilityCone = GEOMETRY_FACTORY.createPolygon(polyCoords.toArray(new Coordinate[0]));
//		Polygon visibilityCone = (Polygon) poly.difference(UnaryUnionOp.union(obstacles));
		visibilityPolygonsCache.put(pair, visibilityCone);
		return visibilityCone;
	}
}
