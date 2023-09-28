/* 
 * Copyright 2023 by Gabriele Filomena
 * University of Liverpool, UK
 * The MIT License (MIT)
 *
 */

package sim.graph;

import java.util.ArrayList;
import java.util.List;
import java.util.Random;

import sim.field.geo.VectorLayer;
import sim.util.geo.MasonGeometry;

/**
 * A class containing functions to identify random nodes, given certain
 * conditions. Usually these are used to identify origin and destination nodes
 * for possible trips.
 *
 */
public class NodesLookup {

    /**
     * Returns a randomly selected node from the given graph.
     *
     * @param network The input graph.
     * @return A randomly selected node from the graph.
     */
    public static NodeGraph randomNode(Graph network) {
		final Random random = new Random();
		final Integer randomInt = random.nextInt(network.nodesMap.values().size());
		final ArrayList<NodeGraph> nodes = new ArrayList<>(network.nodesMap.values());
		return nodes.get(randomInt);
	}

    /**
     * Returns a randomly selected node from the graph that is also contained in the list of node geometries.
     *
     * @param network         The input graph.
     * @param nodesGeometries The list of node geometries.
     * @return A randomly selected node from the graph that meets the specified criteria.
     */
	public static NodeGraph randomNodeFromList(Graph network, ArrayList<MasonGeometry> nodesGeometries) {
		final Random random = new Random();
		final Integer randomInt = random.nextInt(nodesGeometries.size());
		final MasonGeometry geoNode = nodesGeometries.get(randomInt);
		return network.findNode(geoNode.geometry.getCoordinate());
	}

    /**
     * Returns a randomly selected node within a specified radius from the given origin node and outside its region.
     *
     * @param network    The input graph.
     * @param originNode The origin node.
     * @param radius     The maximum distance from the origin node within which the random node should be found.
     * @return A randomly selected node that satisfies the specified conditions.
     */
	public static NodeGraph randomNodeRegion(Graph network, NodeGraph originNode, double radius) {

		final MasonGeometry originNodeGeometry = originNode.masonGeometry;
		Random random = new Random();
		double EXPANSION_FACTOR = 1.10;
		double RADIUS_THRESHOLD = 2.00;
		NodeGraph node = null;
		double expandingRadius = radius;

		while (true) {
			if (expandingRadius >= radius * RADIUS_THRESHOLD)
				return null;

			ArrayList<MasonGeometry> spatialFilter = network.junctions
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
			node = network.findNode(nodeGeometry.geometry.getCoordinate());

			if (node != null) {
				return node;
			}

			expandingRadius *= EXPANSION_FACTOR;
		}
	}

    /**
     * Returns a randomly selected node within a specified radius from the given origin node and outside its region.
     *
     * @param network    The input graph.
     * @param originNode The origin node.
     * @param radius     The maximum distance from the origin node within which the random node should be found.
     * @return A randomly selected node that satisfies the specified conditions.
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
     * Returns a randomly selected node whose distance from the origin node falls within the specified range and has centrality values above or equal to a specified percentile.
     *
     * @param network    The input graph.
     * @param junctions  The vector layer representing junctions.
     * @param originNode The origin node.
     * @param lowerLimit The minimum distance from the origin node.
     * @param upperLimit The maximum distance from the origin node.
     * @param percentile The percentile used as a threshold for centrality values.
     * @return A randomly selected node that satisfies the specified distance and centrality criteria.
     */
	public static ArrayList<NodeGraph> getNodesBetweenDistanceInterval(Graph network, VectorLayer junctions,
			NodeGraph node, double lowerLimit, double upperLimit) {

		ArrayList<NodeGraph> containedNodes = new ArrayList<>();
		MasonGeometry originGeometry = node.masonGeometry;
		ArrayList<MasonGeometry> containedGeometries = junctions.featuresBetweenLimits(originGeometry.geometry,
				lowerLimit, upperLimit);
		for (final MasonGeometry masonGeometry : containedGeometries)
			containedNodes.add(network.findNode(masonGeometry.geometry.getCoordinate()));
		return containedNodes;
	}

    /**
     * Returns a randomly selected node whose distance from the origin node falls within the specified range and has 
     * centrality values above or equal to a specified percentile.
     *
     * @param network    The input graph.
     * @param junctions  The vector layer representing junctions.
     * @param originNode The origin node.
     * @param lowerLimit The minimum distance from the origin node.
     * @param upperLimit The maximum distance from the origin node.
     * @param percentile The percentile used as a threshold for centrality values.
     * @return A randomly selected node that satisfies the specified distance and centrality criteria.
     */
	public static NodeGraph randomNodeBetweenDistanceInterval(Graph network, VectorLayer junctions,
			NodeGraph originNode, double lowerLimit, double upperLimit) {

		final Random random = new Random();
		final ArrayList<NodeGraph> candidates = getNodesBetweenDistanceInterval(network, junctions, originNode,
				lowerLimit, upperLimit);
		final int randomInt = random.nextInt(candidates.size());
		final NodeGraph node = candidates.get(randomInt);
		return node;
	}

	 /**
     * Returns a list of nodes whose distance from the given node falls within the specified range and belong to a different region.
     *
     * @param network    The input graph.
     * @param junctions  The vector layer representing junctions.
     * @param node       The reference node.
     * @param lowerLimit The minimum distance from the reference node.
     * @param upperLimit The maximum distance from the reference node.
     * @return A list of nodes that meet the distance and region criteria.
     */
	public static ArrayList<NodeGraph> getNodesBetweenDistanceIntervalRegion(Graph network, VectorLayer junctions,
			NodeGraph node, double lowerLimit, double upperLimit) {
		ArrayList<NodeGraph> containedNodes = new ArrayList<>();
		containedNodes = getNodesBetweenDistanceInterval(network, junctions, node, lowerLimit, upperLimit);
		return network.nodesInRegion(containedNodes, node.regionID);
	}

	/**
	 * Given a graph, the function returns a random node whose distance from a
	 * passed origin node is within certain limits. The returned node belongs to a
	 * region different from the origin node's region.
	 *
	 * @param network    a graph;
	 * @param originNode a node;
	 * @param lowerLimit the minimum distance from the origin node;
	 * @param upperLimit the maximum distance from the origin node;
	 */
	public static NodeGraph randomNodeBetweenDistanceIntervalRegion(Graph network, VectorLayer junctions,
			NodeGraph originNode, double lowerLimit, double upperLimit) {

		final Random random = new Random();
		final ArrayList<NodeGraph> candidates = getNodesBetweenDistanceIntervalRegion(network, junctions, originNode,
				lowerLimit, upperLimit);
		final int randomInt = random.nextInt(candidates.size());
		final NodeGraph node = candidates.get(randomInt);
		return node;
	}

	/**
	 * Given a node in the graph, it returns all the nodes whose distance from the
	 * given node is between the given lower and the upper limits. Moreover, the
	 * centrality values of these nodes are higher than the value at the passed
	 * percentile. *
	 * 
	 * @param node;
	 * @param lowerLimit;
	 * @param upperLimit;
	 * @param percentile  the percentile to use as threshold;
	 */
	public static ArrayList<NodeGraph> getSalientNodesBetweenDistanceInterval(Graph network, VectorLayer junctions,
			NodeGraph node, double lowerLimit, double upperLimit, double percentile) {

		final ArrayList<NodeGraph> containedNodes = new ArrayList<>();
		final ArrayList<NodeGraph> containedSalientNodes = new ArrayList<>();
		final MasonGeometry originGeometry = node.masonGeometry;
		final ArrayList<MasonGeometry> containedGeometries = junctions.featuresBetweenLimits(originGeometry.geometry,
				lowerLimit, upperLimit);
		for (final MasonGeometry masonGeometry : containedGeometries)
			containedNodes.add(network.findNode(masonGeometry.geometry.getCoordinate()));
		final ArrayList<NodeGraph> salientNodes = new ArrayList<>(network.graphSalientNodes(percentile).keySet());

		for (final NodeGraph otherNode : containedNodes)
			if (salientNodes.contains(otherNode))
				containedSalientNodes.add(otherNode);
		return containedSalientNodes;
	}

    /**
     * Returns a randomly selected node whose distance from the origin node falls within the specified range and has 
     * centrality values above or equal to a specified percentile.
     *
     * @param network    The input graph.
     * @param junctions  The vector layer representing junctions.
     * @param originNode The origin node.
     * @param lowerLimit The minimum distance from the origin node.
     * @param upperLimit The maximum distance from the origin node.
     * @param percentile The percentile used as a threshold for centrality values.
     * @return A randomly selected node that satisfies the specified distance and centrality criteria.
     */
	public static NodeGraph randomSalientNodeBetweenDistanceInterval(Graph network, VectorLayer junctions,
			NodeGraph originNode, double lowerLimit, double upperLimit, double percentile) {

		double PERCENTILE_DECREASE = 0.05;
		final Random random = new Random();
		NodeGraph node = null;
		while (node == null) {
			final ArrayList<NodeGraph> candidates = getSalientNodesBetweenDistanceInterval(network, junctions,
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
     * Returns a randomly selected node whose distance from the origin node falls within the specified range and has 
     * centrality values above or equal to a specified percentile.
     *
     * @param network    The input graph.
     * @param junctions  The vector layer representing junctions.
     * @param originNode The origin node.
     * @param lowerLimit The minimum distance from the origin node.
     * @param upperLimit The maximum distance from the origin node.
     * @param percentile The percentile used as a threshold for centrality values.
     * @return A randomly selected node that satisfies the specified distance and centrality criteria.
     */
	public static NodeGraph randomNodeBetweenDistanceIntervalDMA(Graph network, VectorLayer junctions,
			NodeGraph originNode, double lowerLimit, double upperLimit, String DMA) {

		final Random random = new Random();
		int randomInt;
		NodeGraph node = null;
		double TOLERANCE = 50.0;
		double DISTANCE_MULTIPLIER = 1.50;
		ArrayList<NodeGraph> candidates = new ArrayList<>();

		while (node == null) {
			candidates = getNodesBetweenDistanceInterval(network, junctions, originNode, lowerLimit, upperLimit);
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
     * based on the provided DMA (District Market Area).
     *
     * @param network    The input graph.
     * @param DMA        The desired node category.
     * @return A randomly selected node from the specified category based on the DMA.
     */
	public static NodeGraph randomNodeDMA(Graph network, String DMA) {

		final Random random = new Random();
		ArrayList<NodeGraph> candidates = network.getNodes();

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
