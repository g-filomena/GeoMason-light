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
 * A class containing functions to identify random nodes, given certain
 * conditions. These are used to identify origin and destination nodes
 * for possible trips.
 *
 */
public class NodesLookup {

    /**
     * Returns a randomly selected node from the given graph.
     *
     * @param graph The input graph.
     * @return A randomly selected node from the graph.
     */
    public static NodeGraph randomNode(Graph graph) {
		final Random random = new Random();
		final Integer randomInt = random.nextInt(graph.nodesMap.values().size());
		final ArrayList<NodeGraph> nodes = new ArrayList<>(graph.nodesMap.values());
		return nodes.get(randomInt);
	}

    /**
     * Returns a randomly selected node from the graph that is also contained in a list of node geometries.
     *
     * @param graph           The input graph.
     * @param nodesGeometries The list of node geometries.
     * @return A randomly selected node from the graph that meets the specified criteria.
     */
	public static NodeGraph randomNodeFromList(Graph graph, ArrayList<MasonGeometry> nodesGeometries) {
		final Random random = new Random();
		final Integer randomInt = random.nextInt(nodesGeometries.size());
		final MasonGeometry geoNode = nodesGeometries.get(randomInt);
		return graph.findNode(geoNode.geometry.getCoordinate());
	}

    /**
     * Returns a randomly selected node within a specified radius from the given origin node and outside its region.
     *
     * @param graph      The input graph.
     * @param originNode The origin node.
     * @param radius     The maximum distance from the origin node within which the random node should be found.
     * @return A randomly selected node that satisfies the specified conditions.
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

			ArrayList<MasonGeometry> spatialFilter = graph.junctions
					.featuresWithinDistance(originNodeGeometry.geometry, expandingRadius);
			if (spatialFilter.size() < 1) {
				expandingRadius *= EXPANSION_FACTOR;
				continue;
			}
			VectorLayer junctionsWithin = new VectorLayer(spatialFilter);

			if (junctionsWithin.geometriesList.isEmpty()) {
				expandingRadius *= EXPANSION_FACTOR;
				continue;
			}

			ArrayList<MasonGeometry> regionFilter = junctionsWithin.filterFeatures("district", originNode.regionID,
					false);
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
	 * Returns a randomly selected node from the graph whose distance from the origin node matches one of the provided 
	 * distances.
	 *
	 * @param graph      The input graph.
	 * @param junctions  The vector layer representing junctions.
	 * @param originNode The origin node.
	 * @param distances  The list of possible distances used to identify the node.
	 * @return A randomly selected node that matches the specified distance criteria.
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
		ArrayList<NodeGraph> candidates = new ArrayList<>();

		while (true) {
			final double lowerLimit = distance - TOLERANCE;
			final double upperLimit = distance + TOLERANCE;
			candidates = getNodesBetweenDistanceInterval(graph, junctions, originNode, lowerLimit, upperLimit);
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
	 * Returns a list of nodes in the graph whose distance from the given node falls within the specified range.
	 *
	 * @param graph      The input graph.
	 * @param junctions  The vector layer representing junctions.
	 * @param node       The reference node.
	 * @param lowerLimit The minimum distance from the reference node.
	 * @param upperLimit The maximum distance from the reference node.
	 * @return A list of nodes that meet the distance criteria.
	 */
	public static ArrayList<NodeGraph> getNodesBetweenDistanceInterval(Graph graph, VectorLayer junctions,
			NodeGraph node, double lowerLimit, double upperLimit) {

		ArrayList<NodeGraph> containedNodes = new ArrayList<>();
		MasonGeometry originGeometry = node.masonGeometry;
		ArrayList<MasonGeometry> containedGeometries = junctions.featuresBetweenLimits(originGeometry.geometry,
				lowerLimit, upperLimit);
		for (final MasonGeometry masonGeometry : containedGeometries)
			containedNodes.add(graph.findNode(masonGeometry.geometry.getCoordinate()));
		return containedNodes;
	}

    /**
     * Returns a randomly selected node whose distance from the origin node falls within the specified range.
     *
     * @param graph      The input graph.
     * @param junctions  The vector layer representing junctions.
     * @param originNode The origin node.
     * @param lowerLimit The minimum distance from the origin node.
     * @param upperLimit The maximum distance from the origin node.
     * @return A randomly selected node that satisfies the specified distance and centrality criteria.
     */
	public static NodeGraph randomNodeBetweenDistanceInterval(Graph graph, VectorLayer junctions,
			NodeGraph originNode, double lowerLimit, double upperLimit) {

		final Random random = new Random();
		final ArrayList<NodeGraph> candidates = getNodesBetweenDistanceInterval(graph, junctions, originNode,
				lowerLimit, upperLimit);
		final int randomInt = random.nextInt(candidates.size());
		final NodeGraph node = candidates.get(randomInt);
		return node;
	}

	 /**
     * Returns a list of nodes whose distance from the given node falls within the specified range and belong to a 
     * different region.
     *
     * @param graph      The input graph.
     * @param junctions  The vector layer representing junctions.
     * @param node       The reference node.
     * @param lowerLimit The minimum distance from the reference node.
     * @param upperLimit The maximum distance from the reference node.
     * @return A list of nodes that meet the distance and region criteria.
     */
	public static ArrayList<NodeGraph> getNodesBetweenDistanceIntervalRegion(Graph graph, VectorLayer junctions,
			NodeGraph node, double lowerLimit, double upperLimit) {
		ArrayList<NodeGraph> containedNodes = new ArrayList<>();
		containedNodes = getNodesBetweenDistanceInterval(graph, junctions, node, lowerLimit, upperLimit);
		return graph.nodesInRegion(containedNodes, node.regionID);
	}

	/**
	 * Returns a randomly selected node from the graph whose distance from the origin node falls within the specified
	 * range and belongs to a different region from the origin node's region.
	 *
	 * @param graph      The input graph.
	 * @param junctions  The vector layer representing junctions.
	 * @param originNode The origin node.
	 * @param lowerLimit The minimum distance from the origin node.
	 * @param upperLimit The maximum distance from the origin node.
	 * @return A randomly selected node that satisfies the specified distance and region criteria.
	 */
	public static NodeGraph randomNodeBetweenDistanceIntervalRegion(Graph graph, VectorLayer junctions,
			NodeGraph originNode, double lowerLimit, double upperLimit) {

		final Random random = new Random();
		final ArrayList<NodeGraph> candidates = getNodesBetweenDistanceIntervalRegion(graph, junctions, originNode,
				lowerLimit, upperLimit);
		final int randomInt = random.nextInt(candidates.size());
		final NodeGraph node = candidates.get(randomInt);
		return node;
	}

	/**
	 * Returns a list of nodes in the graph whose distance from the given node falls within the specified range and have
	 * centrality values above the specified percentile.
	 *
	 * @param graph      The input graph.
	 * @param junctions  The vector layer representing junctions.
	 * @param node       The reference node.
	 * @param lowerLimit The minimum distance from the reference node.
	 * @param upperLimit The maximum distance from the reference node.
	 * @param percentile The percentile used as a threshold for centrality values.
	 * @return A list of nodes that meet the distance and centrality criteria.
	 */
	public static ArrayList<NodeGraph> getSalientNodesBetweenDistanceInterval(Graph graph, VectorLayer junctions,
			NodeGraph node, double lowerLimit, double upperLimit, double percentile) {

		final ArrayList<NodeGraph> containedNodes = new ArrayList<>();
		final ArrayList<NodeGraph> containedSalientNodes = new ArrayList<>();
		final MasonGeometry originGeometry = node.masonGeometry;
		final ArrayList<MasonGeometry> containedGeometries = junctions.featuresBetweenLimits(originGeometry.geometry,
				lowerLimit, upperLimit);
		for (final MasonGeometry masonGeometry : containedGeometries)
			containedNodes.add(graph.findNode(masonGeometry.geometry.getCoordinate()));
		final ArrayList<NodeGraph> salientNodes = new ArrayList<>(graph.graphSalientNodes(percentile).keySet());

		for (final NodeGraph otherNode : containedNodes)
			if (salientNodes.contains(otherNode))
				containedSalientNodes.add(otherNode);
		return containedSalientNodes;
	}

    /**
     * Returns a randomly selected node whose distance from the origin node falls within the specified range and has 
     * centrality values above or equal to a specified percentile.
     *
     * @param graph      The input graph.
     * @param junctions  The vector layer representing junctions.
     * @param originNode The origin node.
     * @param lowerLimit The minimum distance from the origin node.
     * @param upperLimit The maximum distance from the origin node.
     * @param percentile The percentile used as a threshold for centrality values.
     * @return A randomly selected node that satisfies the specified distance and centrality criteria.
     */
	public static NodeGraph randomSalientNodeBetweenDistanceInterval(Graph graph, VectorLayer junctions,
			NodeGraph originNode, double lowerLimit, double upperLimit, double percentile) {

		double PERCENTILE_DECREASE = 0.05;
		final Random random = new Random();
		NodeGraph node = null;
		while (node == null) {
			final ArrayList<NodeGraph> candidates = getSalientNodesBetweenDistanceInterval(graph, junctions,
					originNode, lowerLimit, upperLimit, percentile);
			final int randomInt = random.nextInt(candidates.size());
			node = candidates.get(randomInt);
			percentile -= PERCENTILE_DECREASE;
			if (percentile == 0.0)
				return null;
		}
		return node;
	}

	/**
	 * Returns a randomly selected node from the graph whose distance from the origin node falls within the specified 
	 * range and belongs to a certain category ("live," "work," "visit") based on the provided DMA label.
	 * 
	 * @param graph      The input graph.
	 * @param junctions  The vector layer representing junctions.
	 * @param originNode The origin node.
	 * @param lowerLimit The minimum distance from the origin node.
	 * @param upperLimit The maximum distance from the origin node.
	 * @param DMA        The desired node category ("live," "work," "visit") or "random".
	 * @return A randomly selected node that satisfies the specified distance and category criteria.
	 */
	public static NodeGraph randomNodeBetweenDistanceIntervalDMA(Graph graph, VectorLayer junctions,
			NodeGraph originNode, double lowerLimit, double upperLimit, String DMA) {

		final Random random = new Random();
		int randomInt;
		NodeGraph node = null;
		double TOLERANCE = 50.0;
		double DISTANCE_MULTIPLIER = 1.50;
		ArrayList<NodeGraph> candidates = new ArrayList<>();

		while (node == null) {
			candidates = getNodesBetweenDistanceInterval(graph, junctions, originNode, lowerLimit, upperLimit);
			ArrayList<NodeGraph> candidatesDMA = new ArrayList<>();
			if (DMA.equals("random"))
				candidatesDMA = new ArrayList<>(candidates);
			else
				for (final NodeGraph otherNode : candidates)
					if (otherNode.DMA.equals(DMA))
						candidatesDMA.add(otherNode);

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
     * Returns a randomly selected node from the graph that belongs to a specified category ("live," "work," "visit") 
     * based on the provided DMA label.
     *
     * @param graph      The input graph.
	 * @param DMA        The desired node category ("live," "work," "visit") or "random".
     * @return A randomly selected node from the specified category based on the DMA.
     */
	public static NodeGraph randomNodeDMA(Graph graph, String DMA) {

		final Random random = new Random();
		ArrayList<NodeGraph> candidates = graph.getNodes();

		ArrayList<NodeGraph> candidatesDMA = new ArrayList<>();
		if (DMA.equals("random"))
			candidatesDMA = new ArrayList<>(candidates);
		else
			for (final NodeGraph n : candidates)
				if (n.DMA.equals(DMA))
					candidatesDMA.add(n);

		final int randomInt = random.nextInt(candidatesDMA.size());
		return candidatesDMA.get(randomInt);
	}
}
