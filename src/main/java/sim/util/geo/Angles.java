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
package sim.util.geo;

import org.locationtech.jts.geom.Coordinate;
import org.locationtech.jts.geom.Geometry;
import org.locationtech.jts.geom.GeometryFactory;
import org.locationtech.jts.geom.MultiPoint;

import sim.graph.NodeGraph;
import sim.graph.NodeGraphUtils;

/**
 * The class provides utility methods for working with angles and coordinates.
 * It contains static methods to calculate angles between points, compute angular differences,
 * and perform various angle-related calculations.
 */
public class Angles {

	/**
	 * Computes the dot product between two vectors represented as arrays of doubles.
	 *
	 * @param vectorA The first vector as an array of doubles.
	 * @param vectorB The second vector as an array of doubles.
	 * @return The dot product of the two input vectors.
	 */
	private static double dot(double[] vectorA, double[] vectorB) {
	    return vectorA[0] * vectorB[0] + vectorA[1] * vectorB[1];
	}

	/**
	 * Computes the angle in degrees between two NodeGraph instances.
	 *
	 * This method calculates the angle in degrees between two NodeGraph instances, where
	 * 'originNode' represents the starting point, and 'destinationNode' represents the
	 * ending point. The angle is measured with respect to the positive y-axis.
	 *
	 * @param originNode      The starting NodeGraph.
	 * @param destinationNode The ending NodeGraph.
	 * @return The angle in degrees between the two NodeGraph instances with respect to the y-axis.
	 */
	public static double angle(NodeGraph originNode, NodeGraph destinationNode) {
		final Coordinate origin = originNode.getCoordinate();
		final Coordinate destination = destinationNode.getCoordinate();
		final double[] vectorA = {origin.x - origin.x, origin.y - (origin.y + 2000)};
		final double[] vectorB = {origin.x - destination.x, origin.y - destination.y};
		final double dot_prod = dot(vectorA, vectorB);
		final double magA = Math.pow(dot(vectorA, vectorA), 0.5);
		final double magB = Math.pow(dot(vectorB, vectorB), 0.5);

		final double anglRad = Math.acos(dot_prod / magB / magA);
		double angleDeg = Math.toDegrees(anglRad) % 360;
		if (destination.x < origin.x)
			angleDeg = 180 + (180 - angleDeg);
		return angleDeg;
	}

	/**
	 * Computes the angle in degrees formed by two sets of coordinates with respect to the y-axis.
	 *
	 * This method calculates the angle between two sets of coordinates, where the first set
	 * represents the origin and the second set represents the destination. The angle is measured
	 * with respect to the positive y-axis and returned in degrees.
	 *
	 * @param origin      The coordinates of the origin point.
	 * @param destination The coordinates of the destination point.
	 * @return The angle in degrees between the origin and destination with respect to the y-axis.
	 */
	public static double angle(Coordinate origin, Coordinate destination) {
		final double[] vectorA = {origin.x - origin.x, origin.y - (origin.y + 2000)};
		final double[] vectorB = {origin.x - destination.x, origin.y - destination.y};
		final double dot_prod = dot(vectorA, vectorB);
		final double magA = Math.pow(dot(vectorA, vectorA), 0.5);
		final double magB = Math.pow(dot(vectorB, vectorB), 0.5);

		final double anglRad = Math.acos(dot_prod / magB / magA);
		double angleDeg = Math.toDegrees(anglRad) % 360;
		if (destination.x < origin.x)
			angleDeg = 180 + (180 - angleDeg);
		return angleDeg;
	}

	/**
	 * Calculates the angular difference between two angles, considering their orientations.
	 *
	 * This method computes the absolute angular difference between two angles, taking
	 * into account their orientations in degrees. The result is the smallest angular
	 * difference between the two angles, considering both positive and negative
	 * angular directions.
	 *
	 * @param angleA The first angle in degrees.
	 * @param angleB The second angle in degrees.
	 * @return The absolute angular difference between angleA and angleB in degrees.
	 */
	public static double differenceAngles(Double angleA, Double angleB) {
		double difference = 0;
		// check if same quadrant
		if (angleA <= 180 & angleB <= 180 || angleA > 180 & angleB > 180)
			difference = Math.abs(angleA - angleB);

		else if (angleA > 180 & angleB <= 180) {
			final double tmpA = Math.abs(angleB - angleA);
			final double tmpB = Math.abs(angleA - (angleB + 360));
			difference = Math.min(tmpA, tmpB);
		}
		// (angleB > 180 & angleA <= 180)
		else {
			final double tmpA = Math.abs(angleB - angleA);
			final double tmpB = Math.abs(angleB - (angleA + 360));
			difference = Math.min(tmpA, tmpB);
		}
		return difference;
	}

	/**
	 * Checks whether an angle between an origin and a second node falls within a specified cone.
	 *
	 * This method determines if an angle, given the angle between the origin and
	 * the destination node (angleOD), the angle between the origin and a possible
	 * intermediate node (angleON), and the amplitude of the cone (cone), lies within
	 * the specified cone.
	 *
	 * @param angleOD The angle between the origin and the destination node in degrees.
	 * @param angleON The angle between the origin and the intermediate node in degrees.
	 * @param cone    The amplitude of the cone in degrees.
	 * @return        {@code true} if the angle between the origin and the intermediate
	 *                node is within the specified cone; {@code false} otherwise.
	 */
	public static boolean isInDirection(double angleOD, double angleON, double cone) {
		double limitLeft = angleOD - (cone / 2.0 + 1.0);
		double limitRight = angleOD + (cone / 2.0 + 1.0);

		if (limitLeft < 0.0)
			limitLeft = 360.0 + limitLeft;
		if (limitRight > 360.0)
			limitRight = limitRight - 360.0;

		// over the 0
		if (limitLeft > limitRight) {
			if (angleON >= limitLeft)
				return true;
		} else if (limitLeft < limitRight)
			if (angleON >= limitLeft && angleON <= limitRight)
				return true;
		return false;
	}

	/**
	 * Calculates the coordinate at a specified angle and distance from the given origin node.
	 *
	 * Given a NodeGraph, a distance from it, and a desired angle formed with the y
	 * axis, it returns the coordinates of the resulting location.
	 *
	 * @param originNode The origin node from which to measure the angle and distance.
	 * @param distance   The distance from the origin node to the target coordinate.
	 * @param angle      The angle in degrees at which to position the target coordinate.
	 * @return           The calculated coordinate based on the given parameters.
	 */
	public static Coordinate getCoordAngle(NodeGraph originNode, double distance, double angle) {
		final double x = distance * Math.sin(Math.toRadians(angle)) + originNode.getCoordinate().x;
		final double y = distance * Math.cos(Math.toRadians(angle)) + originNode.getCoordinate().y;
		final Coordinate coord = new Coordinate(x, y);
		return coord;
	}

	/**
	 * Calculates the coordinate of a point within the view cone of a node.
	 *
	 * Given two NodeGraph instances and half of the amplitude of a desired angle (either
	 * positive or negative, e.g. + 35.0° or - 35.0° --> 70°), this method computes the coordinate of a point that, 
	 * when connected with the first node, forms one of the lines of the view cone.
	 *
	 * @param node          The starting node.
	 * @param otherNode     The target node.
	 * @param desiredAngle  The half of the field of view's amplitude in degrees.
	 * @return              The coordinate of the point within the view cone.
	 */
	private static Coordinate angleViewField(NodeGraph node, NodeGraph otherNode, double desiredAngle) {
		double angle = angle(node, otherNode);

		if (angle > 360)
			angle = angle - 360.0;
		double resultingAngle = angle + desiredAngle;
		if (resultingAngle > 360.0)
			resultingAngle = resultingAngle - 360.0;
		final Coordinate coord = getCoordAngle(node, NodeGraphUtils.nodesDistance(node, otherNode), resultingAngle);
		return coord;
	}

	/**
	 * Computes the view field geometry between two nodes based on a specified field of view angle.
	 * The view field represents the region visible from the first node to the second node within
	 * the given field of view angle.
	 *
	 * @param node         The first node, from which the view field is observed.
	 * @param otherNode    The second node, defining the direction of observation.
	 * @param fieldOfView  The field of view angle in degrees (half of the total angle).
	 * @return A geometry representing the view field between the two nodes.
	 */
	public static Geometry viewField(NodeGraph node, NodeGraph otherNode, double fieldOfView) {

		if (fieldOfView >= 180.0)
			fieldOfView = 140.0;

		final Coordinate coordA = angleViewField(node, otherNode, fieldOfView / 2.0);
		final Coordinate coordB = angleViewField(node, otherNode, -fieldOfView / 2.0);
		final GeometryFactory geometryFactory = new GeometryFactory();
		Coordinate[] coordinates = { node.getCoordinate(), coordA, coordB, otherNode.getCoordinate() };
		MultiPoint points = geometryFactory.createMultiPointFromCoords(coordinates);

		final Geometry viewField = points.convexHull();
		return viewField;
	}
}
