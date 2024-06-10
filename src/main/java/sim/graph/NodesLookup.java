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
import java.util.List;
import java.util.Random;

import sim.field.geo.VectorLayer;
import sim.util.geo.MasonGeometry;

/**
 * A class containing functions to identify random nodes within a graph, based
 * on specific criteria such as location, region, and distance. These methods
 * are useful for generating origin and destination nodes for agent-based
 * simulations in spatial models.
 */
public class NodesLookup {

	/**
	 * Returns a randomly selected node from a graph. This method is useful for
	 * generating random starting or ending points in simulations.
	 *
	 * @param graph The graph from which a random node is to be selected.
	 * @return A randomly chosen NodeGraph object from the graph.
	 */
	public static NodeGraph randomNode(Graph graph) {
		final Random random = new Random();
		final Integer randomInt = random.nextInt(graph.nodesGraph.size());
		final List<NodeGraph> nodes = new ArrayList<>(graph.nodesGraph);
		return nodes.get(randomInt);
	}

	/**
	 * Returns a randomly selected node from a specified list of node geometries
	 * within a graph. This method allows for the random selection of nodes that
	 * meet certain spatial criteria represented by the node geometries.
	 *
	 * @param graph           The graph containing the nodes.
	 * @param nodesGeometries A List of MasonGeometry objects representing specific
	 *                        node locations.
	 * @return A NodeGraph object randomly selected from the specified geometries.
	 */
	public static NodeGraph randomNodeFromList(Graph graph, List<MasonGeometry> nodesGeometries) {
		final Random random = new Random();
		final Integer randomInt = random.nextInt(nodesGeometries.size());
		final MasonGeometry geoNode = nodesGeometries.get(randomInt);
		return graph.findNode(geoNode.geometry.getCoordinate());
	}

	/**
	 * Returns a randomly selected node from a graph within a specified radius from
	 * an origin node and outside the origin node's region. This method is used to
	 * generate nodes that are geographically spread yet regionally distinct from
	 * the origin.
	 *
	 * @param graph      The graph from which to select the node.
	 * @param originNode The origin node serving as the center of the search radius.
	 * @param radius     The radius within which to search for a suitable node.
	 * @return A randomly selected NodeGraph object that meets the specified
	 *         criteria.
	 */
	public static NodeGraph randomNodeRegion(Graph graph, NodeGraph originNode, double radius) {

		final MasonGeometry originNodeGeometry = originNode.masonGeometry;
		Random random = new Random();
		double EXPANSION_FACTOR = 1.10;
		double RADIUS_THRESHOLD = 2.00;
		NodeGraph node = null;
		double expandingRadius = radius;

		while (true) {
			if (expandingRadius >= radius * RADIUS_THRESHOLD)
				return null;

			List<MasonGeometry> spatialFilter = graph.junctions.featuresWithinDistance(originNodeGeometry.geometry,
					expandingRadius);
			if (spatialFilter.size() < 1) {
				expandingRadius *= EXPANSION_FACTOR;
				continue;
			}
			VectorLayer junctionsWithin = new VectorLayer(spatialFilter);

			if (junctionsWithin.getGeometries().isEmpty()) {
				expandingRadius *= EXPANSION_FACTOR;
				continue;
			}

			List<MasonGeometry> regionFilter = junctionsWithin.filterFeatures("district", originNode.regionID, false);
			if (regionFilter.isEmpty()) {
				expandingRadius *= EXPANSION_FACTOR;
				continue;
			}

			int randomInt = random.nextInt(regionFilter.size());
			MasonGeometry nodeGeometry = regionFilter.get(randomInt);
			node = graph.findNode(nodeGeometry.geometry.getCoordinate());

			if (node != null)
				return node;

			expandingRadius *= EXPANSION_FACTOR;
		}
	}

	/**
	 * Returns a randomly selected node from a graph whose distance from the origin
	 * node matches one of the provided distances.
	 *
	 * @param graph      The input graph.
	 * @param junctions  The vector layer representing junctions.
	 * @param originNode The origin node.
	 * @param distances  The list of possible distances used to identify the node.
	 * @return A randomly selected node that matches the specified distance
	 *         criteria.
	 */
	public static NodeGraph randomNodeFromDistancesSet(Graph graph, VectorLayer junctions, NodeGraph originNode,
			List<Float> distances) {

		double MIN_DISTANCE = 100;
		double TOLERANCE = 50;
		Random random = new Random();
		int randomInt = random.nextInt(distances.size());
		double distance = distances.get(randomInt);
		if (distance < MIN_DISTANCE)
			distance = MIN_DISTANCE;
		NodeGraph node = null;
		List<NodeGraph> candidates = new ArrayList<>();

		while (true) {
			final double lowerLimit = distance - TOLERANCE;
			final double upperLimit = distance + TOLERANCE;
			candidates = getNodesBetweenDistanceInterval(graph, originNode, lowerLimit, upperLimit);
			if (!candidates.isEmpty())
				break;
			else
				TOLERANCE += TOLERANCE;
		}

		while (node == null || node.getID() == originNode.getID()) {
			randomInt = random.nextInt(candidates.size());
			node = candidates.get(randomInt);
			if (graph.getEdgeBetween(originNode, node) != null)
				node = null;
		}
		return node;
	}

	/**
	 * Returns a List of nodes in the graph whose distance from the given node falls
	 * within the specified range.
	 *
	 * @param graph      The input graph.
	 * @param node       The reference node.
	 * @param lowerLimit The minimum distance from the reference node.
	 * @param upperLimit The maximum distance from the reference node.
	 * @return A List of nodes that meet the distance criteria.
	 */
	public static List<NodeGraph> getNodesBetweenDistanceInterval(Graph graph, NodeGraph node, double lowerLimit,
			double upperLimit) {

		List<NodeGraph> containedNodes = new ArrayList<>();
		MasonGeometry originGeometry = node.masonGeometry;
		List<MasonGeometry> containedGeometries = graph.junctions.featuresBetweenLimits(originGeometry.geometry,
				lowerLimit, upperLimit);
		for (final MasonGeometry masonGeometry : containedGeometries) {
			NodeGraph potentialNode = graph.findNode(masonGeometry.geometry.getCoordinate());
			if (!node.equals(potentialNode))
				containedNodes.add(potentialNode);
		}
		return containedNodes;
	}

	/**
	 * Returns a randomly selected node from a graph whose distance from an origin
	 * node falls within a specified range. This method allows for the selection of
	 * nodes based on proximity constraints.
	 *
	 * @param graph      The graph to search within.
	 * @param originNode The origin node to measure distances from.
	 * @param lowerLimit The minimum distance from the origin node.
	 * @param upperLimit The maximum distance from the origin node.
	 * @return A NodeGraph object randomly selected within the specified distance
	 *         range.
	 */
	public static NodeGraph randomNodeBetweenDistanceInterval(Graph graph, NodeGraph originNode, double lowerLimit,
			double upperLimit) {

		final Random random = new Random();
		final List<NodeGraph> candidates = getNodesBetweenDistanceInterval(graph, originNode, lowerLimit, upperLimit);
		final int randomInt = random.nextInt(candidates.size());
		final NodeGraph node = candidates.get(randomInt);
		return node;
	}

	/**
	 * Returns a List of nodes whose distance from the given node falls within the
	 * specified range and belong to a different region.
	 *
	 * @param graph      The input graph.
	 * @param node       The reference node.
	 * @param lowerLimit The minimum distance from the reference node.
	 * @param upperLimit The maximum distance from the reference node.
	 * @return A List of nodes that meet the distance and region criteria.
	 */
	public static List<NodeGraph> getNodesBetweenDistanceIntervalRegion(Graph graph, NodeGraph node, double lowerLimit,
			double upperLimit) {
		List<NodeGraph> containedNodes = new ArrayList<>();
		containedNodes = getNodesBetweenDistanceInterval(graph, node, lowerLimit, upperLimit);
		return graph.nodesInRegion(containedNodes, node.regionID);
	}

	/**
	 * Returns a randomly selected node from a graph whose distance from an origin
	 * node falls within a specified range and belongs to a different region than
	 * the origin. This method is useful for simulating regional travel or migration
	 * patterns.
	 *
	 * @param graph      The graph to search within.
	 * @param originNode The origin node to measure distances from.
	 * @param lowerLimit The minimum distance from the origin node.
	 * @param upperLimit The maximum distance from the origin node.
	 * @return A NodeGraph object randomly selected within the specified distance
	 *         range and different region.
	 */
	public static NodeGraph randomNodeBetweenDistanceIntervalRegion(Graph graph, NodeGraph originNode,
			double lowerLimit, double upperLimit) {

		final Random random = new Random();
		final List<NodeGraph> candidates = getNodesBetweenDistanceIntervalRegion(graph, originNode, lowerLimit,
				upperLimit);
		final int randomInt = random.nextInt(candidates.size());
		final NodeGraph node = candidates.get(randomInt);
		return node;
	}

	/**
	 * Returns a randomly selected node whose distance from the origin node falls
	 * within the specified range and has centrality values above or equal to a
	 * specified percentile.
	 *
	 * @param graph      The input graph.
	 * @param originNode The origin node.
	 * @param lowerLimit The minimum distance from the origin node.
	 * @param upperLimit The maximum distance from the origin node.
	 * @param percentile The percentile used as a threshold for centrality values.
	 * @return A randomly selected node that satisfies the specified distance and
	 *         centrality criteria.
	 */
	public static NodeGraph randomSalientNodeBetweenDistanceInterval(Graph graph, NodeGraph originNode,
			double lowerLimit, double upperLimit, double percentile) {

		double PERCENTILE_DECREASE = 0.05;
		final Random random = new Random();
		NodeGraph node = null;
		while (node == null) {
			final List<NodeGraph> candidates = graph.getSalientNodesBetweenDistanceInterval(originNode, lowerLimit,
					upperLimit, percentile);
			final int randomInt = random.nextInt(candidates.size());
			node = candidates.get(randomInt);
			percentile -= PERCENTILE_DECREASE;
			if (percentile == 0.0)
				return null;
		}
		return node;
	}

	/**
	 * Returns a random node from the graph whose distance from the origin node
	 * falls within a specified range and belongs to a specific category based on
	 * the DMA label. This method is useful for selecting nodes that are relevant to
	 * specific activities or functions like living, working, or visiting.
	 *
	 * @param graph      The graph to search within.
	 * @param originNode The origin node to measure distances from.
	 * @param lowerLimit The minimum distance from the origin node.
	 * @param upperLimit The maximum distance from the origin node.
	 * @param DMA        The category label (e.g., "live," "work," "visit") or
	 *                   "random" for any category.
	 * @return A NodeGraph object randomly selected based on the specified distance
	 *         and category criteria.
	 */

	public static NodeGraph randomNodeBetweenDistanceIntervalDMA(Graph graph, NodeGraph originNode, double lowerLimit,
			double upperLimit, String DMA) {

		final Random random = new Random();
		int randomInt;
		NodeGraph node = null;
		double TOLERANCE = 50.0;
		double DISTANCE_MULTIPLIER = 1.50;
		List<NodeGraph> candidates = new ArrayList<>();

		while (node == null) {
			candidates = getNodesBetweenDistanceInterval(graph, originNode, lowerLimit, upperLimit);
			List<NodeGraph> candidatesDMA = new ArrayList<>();
			if (DMA.equals("random"))
				candidatesDMA = new ArrayList<>(candidates);
			else {
				for (final NodeGraph otherNode : candidates)
					if (otherNode.DMA.equals(DMA))
						candidatesDMA.add(otherNode);
			}
			if (candidatesDMA.isEmpty()) {
				if (upperLimit > upperLimit * DISTANCE_MULTIPLIER) {
					randomInt = random.nextInt(candidates.size());
					node = candidates.get(randomInt);
					break;
				}
				upperLimit += TOLERANCE;
				continue;
			}

			randomInt = random.nextInt(candidatesDMA.size());
			node = candidatesDMA.get(randomInt);
		}
		return node;
	}

	/**
	 * Returns a randomly selected node from the graph that belongs to a specified
	 * category ("live," "work," "visit") based on the provided DMA label.
	 *
	 * @param graph The input graph.
	 * @param DMA   The desired node category ("live," "work," "visit") or "random".
	 * @return A randomly selected node from the specified category based on the
	 *         DMA.
	 */
	public static NodeGraph randomNodeDMA(Graph graph, String DMA) {

		final Random random = new Random();
		List<NodeGraph> candidates = graph.getNodes();

		List<NodeGraph> candidatesDMA = new ArrayList<>();
		if (DMA.equals("random"))
			candidatesDMA = new ArrayList<>(candidates);
		else
			for (final NodeGraph node : candidates)
				if (node.DMA.equals(DMA))
					candidatesDMA.add(node);

		final int randomInt = random.nextInt(candidatesDMA.size());
		return candidatesDMA.get(randomInt);
	}
}