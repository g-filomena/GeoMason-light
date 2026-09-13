/*
 * Copyright (c) 2023 Gabriele Filomena University of Liverpool, UK
 *
 * This program is free software: it can redistributed and/or modified under the terms of the GNU
 * General Public License 3.0 as published by the Free Software Foundation.
 *
 * See the file "LICENSE" for more information
 */
package sim.graph;

import java.util.ArrayList;
import java.util.List;
import ec.util.MersenneTwisterFast;
import java.util.stream.Collectors;

import sim.field.geo.VectorLayer;

import sim.util.geo.MasonGeometry;

/**
 * A class containing functions to identify random nodes within a graph, based
 * on specific criteria such as location, region, and distance. These methods
 * are useful for generating origin and destination nodes for agent-based
 * simulations in spatial models.
 *
 * <p>
 * Selection methods return {@code null} when no candidate satisfies the
 * criteria; they never throw on empty candidate sets.
 *
 * <p>
 * Every drawing method comes in two forms. The one without a generator draws
 * from a per-thread {@link MersenneTwisterFast}: contention-free under
 * parallel simulation, and not reproducible, which is what these lookups have
 * always been. The overload taking a {@link MersenneTwisterFast} draws from it instead, so a caller
 * that owns a seeded generator - one per agent, say, seeded from the model's
 * seed - gets the same nodes on every run of the same seed, concurrency
 * included. A simulation whose origins and destinations come from the first
 * form cannot be replayed, however carefully everything else is seeded.
 */
public class NodesLookup {

	/**
	 * The generator used by the methods that are not given one: one per thread, so nothing is
	 * shared and nothing contends, seeded from the clock and therefore not reproducible. That is
	 * what these lookups have always been; pass a generator to get repeatability.
	 */
	private static final ThreadLocal<MersenneTwisterFast> FALLBACK =
			ThreadLocal.withInitial(MersenneTwisterFast::new);

	private static MersenneTwisterFast fallbackGenerator() {
		return FALLBACK.get();
	}

	final static double PERCENTILE_DECREASE = 0.05;
	final static double EXPANSION_FACTOR = 1.10;
	final static double RADIUS_THRESHOLD = 2.00;

	final static double MIN_DISTANCE = 100;
	final static double INITIAL_TOLERANCE = 50;
	final static double TOLERANCE_INCREMENT = 50;
	final static double DISTANCE_MULTIPLIER = 1.50;

	/** Hard cap on search-interval expansions before a lookup gives up. */
	final static int MAX_EXPANSIONS = 1000;

	/** Hard cap on random re-draws before a lookup gives up. */
	final static int MAX_DRAW_ATTEMPTS = 100;

	/**
	 * Returns a randomly selected node from a graph. This method is useful for
	 * generating random starting or ending points in simulations.
	 *
	 * @param graph The graph from which a random node is to be selected.
	 * @return A randomly chosen NodeGraph object from the graph.
	 */
	public static NodeGraph randomNode(Graph graph) {
		return randomNode(graph, fallbackGenerator());
	}

	/**
	 * As {@link #randomNode(Graph)}, drawing from the supplied generator.
	 *
	 * @param graph  The graph from which a random node is to be selected.
	 * @param random The generator to draw from.
	 * @return A randomly chosen NodeGraph object from the graph.
	 */
	public static NodeGraph randomNode(Graph graph, MersenneTwisterFast random) {
		final List<NodeGraph> candidates = new ArrayList<>(graph.nodesGraph);
		return selectRandomNode(candidates, random);
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
	public static NodeGraph randomNodeFromGeometries(Graph graph, List<MasonGeometry> nodesGeometries) {
		return randomNodeFromGeometries(graph, nodesGeometries, fallbackGenerator());
	}

	/**
	 * As {@link #randomNodeFromGeometries(Graph, List)}, drawing from the supplied
	 * generator.
	 *
	 * @param graph           The graph containing the nodes.
	 * @param nodesGeometries A List of MasonGeometry objects representing specific
	 *                        node locations.
	 * @param random          The generator to draw from.
	 * @return A NodeGraph object randomly selected from the specified geometries.
	 */
	public static NodeGraph randomNodeFromGeometries(Graph graph, List<MasonGeometry> nodesGeometries,
			MersenneTwisterFast random) {
		if (nodesGeometries.isEmpty()) {
			return null;
		}
		Integer randomInt = random.nextInt(nodesGeometries.size());
		MasonGeometry geoNode = nodesGeometries.get(randomInt);
		return graph.findNode(geoNode.geometry.getCoordinate());
	}

	/**
	 * Returns a randomly selected node lying within a radius of the origin node
	 * and in the same region as it; use
	 * {@link #randomNodeBetweenDistanceIntervalRegion(Graph, NodeGraph, double, double)}
	 * for a node in a different region. The radius grows until a candidate is
	 * found or it has doubled.
	 *
	 * @param graph      The graph from which to select the node.
	 * @param originNode The origin node serving as the center of the search radius.
	 * @param radius     The radius within which to search for a suitable node.
	 * @return A randomly selected NodeGraph object that meets the specified
	 *         criteria.
	 */
	public static NodeGraph randomNodeRegion(Graph graph, NodeGraph originNode, double radius) {
		return randomNodeRegion(graph, originNode, radius, fallbackGenerator());
	}

	/**
	 * As {@link #randomNodeRegion(Graph, NodeGraph, double)}, drawing from the
	 * supplied generator.
	 *
	 * @param graph      The graph from which to select the node.
	 * @param originNode The origin node serving as the center of the search radius.
	 * @param radius     The radius within which to search for a suitable node.
	 * @param random     The generator to draw from.
	 * @return A randomly selected NodeGraph object that meets the specified
	 *         criteria.
	 */
	public static NodeGraph randomNodeRegion(Graph graph, NodeGraph originNode, double radius, MersenneTwisterFast random) {

		final MasonGeometry originNodeGeometry = originNode.masonGeometry;
		double expandingRadius = radius;

		// The region is read off the node, as it is in getNodesBetweenDistanceIntervalRegion.
		// Do not read it from a "district" attribute on the junction geometry: a node generated
		// from a segment endpoint carries no attributes at all.
		while (expandingRadius < radius * RADIUS_THRESHOLD) {

			List<MasonGeometry> spatialFilter = graph.junctions.featuresWithinDistance(originNodeGeometry.geometry,
					expandingRadius);
			List<NodeGraph> candidates = new ArrayList<>();

			for (MasonGeometry junction : spatialFilter) {
				NodeGraph node = graph.findNode(junction.geometry.getCoordinate());
				if (node != null && !node.equals(originNode) && node.getRegionID() == originNode.getRegionID()) {
					candidates.add(node);
				}
			}

			if (!candidates.isEmpty()) {
				return selectRandomNode(candidates, random);
			}
			expandingRadius *= EXPANSION_FACTOR;
		}
		return null;
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
		return randomNodeFromDistancesSet(graph, junctions, originNode, distances, fallbackGenerator());
	}

	/**
	 * As {@link #randomNodeFromDistancesSet(Graph, VectorLayer, NodeGraph, List)},
	 * drawing from the supplied generator.
	 *
	 * @param graph      The input graph.
	 * @param junctions  The vector layer representing junctions.
	 * @param originNode The origin node.
	 * @param distances  The list of possible distances used to identify the node.
	 * @param random     The generator to draw from.
	 * @return A randomly selected node that matches the specified distance
	 *         criteria.
	 */
	public static NodeGraph randomNodeFromDistancesSet(Graph graph, VectorLayer junctions, NodeGraph originNode,
			List<Float> distances, MersenneTwisterFast random) {

		if (distances.isEmpty()) {
			return null;
		}

		// Select a random distance from the list
		double distance = distances.get(random.nextInt(distances.size()));
		if (distance < MIN_DISTANCE) {
			distance = MIN_DISTANCE;
		}

		List<NodeGraph> candidates = new ArrayList<>();
		double tolerance = INITIAL_TOLERANCE;

		for (int expansion = 0; candidates.isEmpty(); expansion++) {
			if (expansion >= MAX_EXPANSIONS) {
				return null;
			}
			double lowerLimit = distance - tolerance;
			double upperLimit = distance + tolerance;
			candidates = getNodesBetweenDistanceInterval(graph, originNode, lowerLimit, upperLimit);
			tolerance += TOLERANCE_INCREMENT;
		}

		// Prefer a candidate that is not directly connected to the origin; bounded so that a
		// candidate set made entirely of neighbours cannot spin this loop forever.
		for (int attempt = 0; attempt < MAX_DRAW_ATTEMPTS; attempt++) {
			NodeGraph node = candidates.get(random.nextInt(candidates.size()));
			if (node.getID() != originNode.getID() && graph.getEdgeBetween(originNode, node) == null) {
				return node;
			}
		}
		return selectRandomNode(candidates, random);
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
		MasonGeometry originGeometry = node.masonGeometry;
		List<MasonGeometry> containedGeometries = graph.junctions.featuresBetweenLimits(originGeometry.geometry,
				lowerLimit, upperLimit);
		return containedGeometries.stream().map(masonGeometry -> graph.findNode(masonGeometry.geometry.getCoordinate()))
				.filter(potentialNode -> !node.equals(potentialNode)).collect(Collectors.toList());
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
		return randomNodeBetweenDistanceInterval(graph, originNode, lowerLimit, upperLimit,
				fallbackGenerator());
	}

	/**
	 * As {@link #randomNodeBetweenDistanceInterval(Graph, NodeGraph, double, double)},
	 * drawing from the supplied generator.
	 *
	 * @param graph      The graph to search within.
	 * @param originNode The origin node to measure distances from.
	 * @param lowerLimit The minimum distance from the origin node.
	 * @param upperLimit The maximum distance from the origin node.
	 * @param random     The generator to draw from.
	 * @return A NodeGraph object randomly selected within the specified distance
	 *         range.
	 */
	public static NodeGraph randomNodeBetweenDistanceInterval(Graph graph, NodeGraph originNode, double lowerLimit,
			double upperLimit, MersenneTwisterFast random) {

		final List<NodeGraph> candidates = getNodesBetweenDistanceInterval(graph, originNode, lowerLimit, upperLimit);
		return selectRandomNode(candidates, random);
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
		return randomNodeBetweenDistanceIntervalRegion(graph, originNode, lowerLimit, upperLimit,
				fallbackGenerator());
	}

	/**
	 * As
	 * {@link #randomNodeBetweenDistanceIntervalRegion(Graph, NodeGraph, double, double)},
	 * drawing from the supplied generator.
	 *
	 * @param graph      The graph to search within.
	 * @param originNode The origin node to measure distances from.
	 * @param lowerLimit The minimum distance from the origin node.
	 * @param upperLimit The maximum distance from the origin node.
	 * @param random     The generator to draw from.
	 * @return A NodeGraph object randomly selected within the specified distance
	 *         range and different region.
	 */
	public static NodeGraph randomNodeBetweenDistanceIntervalRegion(Graph graph, NodeGraph originNode,
			double lowerLimit, double upperLimit, MersenneTwisterFast random) {

		List<NodeGraph> candidates = getNodesBetweenDistanceIntervalRegion(graph, originNode, lowerLimit, upperLimit);
		return selectRandomNode(candidates, random);
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
		return randomSalientNodeBetweenDistanceInterval(graph, originNode, lowerLimit, upperLimit, percentile,
				fallbackGenerator());
	}

	/**
	 * As
	 * {@link #randomSalientNodeBetweenDistanceInterval(Graph, NodeGraph, double, double, double)},
	 * drawing from the supplied generator.
	 *
	 * @param graph      The input graph.
	 * @param originNode The origin node.
	 * @param lowerLimit The minimum distance from the origin node.
	 * @param upperLimit The maximum distance from the origin node.
	 * @param percentile The percentile used as a threshold for centrality values.
	 * @param random     The generator to draw from.
	 * @return A randomly selected node that satisfies the specified distance and
	 *         centrality criteria.
	 */
	public static NodeGraph randomSalientNodeBetweenDistanceInterval(Graph graph, NodeGraph originNode,
			double lowerLimit, double upperLimit, double percentile, MersenneTwisterFast random) {

		// Lower the centrality bar until the interval yields a candidate: an empty candidate set
		// retries at a lower percentile, and a set that yields a node returns it immediately rather
		// than falling through to the next decrement.
		while (percentile > 0.0) {
			List<NodeGraph> candidates = graph.getSalientNodesBetweenDistanceInterval(originNode, lowerLimit,
					upperLimit, percentile);

			if (!candidates.isEmpty()) {
				return selectRandomNode(candidates, random);
			}
			percentile -= PERCENTILE_DECREASE;
		}
		return null;
	}

	/**
	 * Returns a random node from the graph whose distance from the origin node
	 * falls within a specified range and belongs to a specific category based on
	 * the DMA label. This method is useful for selecting nodes that are relevant to
	 * specific activities or functions like living, working, or visiting.
	 *
	 * DMA = Urban Density (D), Mix(M) and Access (A).
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
		return randomNodeBetweenDistanceIntervalDMA(graph, originNode, lowerLimit, upperLimit, DMA,
				fallbackGenerator());
	}

	/**
	 * As
	 * {@link #randomNodeBetweenDistanceIntervalDMA(Graph, NodeGraph, double, double, String)},
	 * drawing from the supplied generator.
	 *
	 * @param graph      The graph to search within.
	 * @param originNode The origin node to measure distances from.
	 * @param lowerLimit The minimum distance from the origin node.
	 * @param upperLimit The maximum distance from the origin node.
	 * @param DMA        The category label, or "random" for any category.
	 * @param random     The generator to draw from.
	 * @return A NodeGraph object randomly selected based on the specified distance
	 *         and category criteria.
	 */
	public static NodeGraph randomNodeBetweenDistanceIntervalDMA(Graph graph, NodeGraph originNode, double lowerLimit,
			double upperLimit, String DMA, MersenneTwisterFast random) {

		// DMA filtering applies until the interval has been widened past this cap; beyond it, any
		// candidate is accepted. The previous cap condition (upperLimit > upperLimit * multiplier)
		// was always false, so a graph without matching DMA nodes span this loop forever.
		final double maxUpperLimitDMA = upperLimit * DISTANCE_MULTIPLIER;

		for (int expansion = 0; expansion < MAX_EXPANSIONS; expansion++) {
			List<NodeGraph> candidates = getNodesBetweenDistanceInterval(graph, originNode, lowerLimit, upperLimit);

			if (upperLimit <= maxUpperLimitDMA) {
				List<NodeGraph> candidatesDMA = getCandidatesByDMA(candidates, DMA);
				if (!candidatesDMA.isEmpty()) {
					return selectRandomNode(candidatesDMA, random);
				}
			} else if (!candidates.isEmpty()) {
				return selectRandomNode(candidates, random);
			}
			// Widen both ends. Moving only the upper one, as this did, meant a graph too sparse to
			// answer the first interval could only ever be answered by a node further away than
			// asked for - a bias in one direction, produced by the search rather than by the data.
			lowerLimit = Math.max(0.0, lowerLimit - INITIAL_TOLERANCE);
			upperLimit += INITIAL_TOLERANCE;
		}
		return null;
	}

	/**
	 * Returns a randomly selected node from the graph that belongs to a specified
	 * category ("live," "work," "visit") based on the provided DMA label.
	 *
	 * DMA = Urban Density (D), Mix(M) and Access (A).
	 *
	 * @param graph The input graph.
	 * @param DMA   The desired node category ("live," "work," "visit") or "random".
	 * @return A randomly selected node from the specified category based on the
	 *         DMA.
	 */
	public static NodeGraph randomNodeDMA(Graph graph, String DMA) {
		return randomNodeDMA(graph, DMA, fallbackGenerator());
	}

	/**
	 * As {@link #randomNodeDMA(Graph, String)}, drawing from the supplied
	 * generator.
	 *
	 * @param graph  The input graph.
	 * @param DMA    The desired node category ("live," "work," "visit") or "random".
	 * @param random The generator to draw from.
	 * @return A randomly selected node from the specified category based on the
	 *         DMA.
	 */
	public static NodeGraph randomNodeDMA(Graph graph, String DMA, MersenneTwisterFast random) {
		List<NodeGraph> candidates = graph.getNodes();
		List<NodeGraph> candidatesDMA = getCandidatesByDMA(candidates, DMA);
		return selectRandomNode(candidatesDMA, random);
	}

	/**
	 * Filters a list of nodes based on the specified DMA. - If DMA is "random",
	 * returns nodes whose DMA is either "work", "visit", or "live". - If DMA is
	 * "workOrVisit", returns nodes whose DMA is either "work" or "visit". -
	 * Otherwise, returns nodes whose DMA exactly matches the given DMA string.
	 *
	 * DMA = Urban Density (D), Mix(M) and Access (A).
	 *
	 * @param nodes List of nodes to filter.
	 * @param DMA   The DMA criterion used for filtering ("random", "workOrVisit",
	 *              or specific DMA).
	 * @return A filtered list of nodes matching the DMA criteria.
	 */
	public static List<NodeGraph> getCandidatesByDMA(List<NodeGraph> nodes, String DMA) {
		if (DMA.equals("random")) {
			return nodes.stream()
					.filter(node -> node.dma.equals("work") || node.dma.equals("visit") || node.dma.equals("live"))
					.collect(Collectors.toList());
		} else if (DMA.equals("workOrVisit")) {
			return nodes.stream().filter(node -> node.dma.equals("work") || node.dma.equals("visit"))
					.collect(Collectors.toList());
		} else {
			return nodes.stream().filter(node -> node.dma.equals(DMA)).collect(Collectors.toList());
		}
	}

	/**
	 * Retrieves nodes from a graph based on the specified DMA.
	 *
	 * DMA = Urban Density (D), Mix(M) and Access (A).
	 *
	 * @param graph The graph containing the nodes to filter.
	 * @param DMA   The DMA criterion used for filtering ("random", "workOrVisit",
	 *              or specific DMA).
	 * @return A filtered list of nodes matching the DMA criteria.
	 */
	public static List<NodeGraph> getNodesByDMA(Graph graph, String DMA) {
		return getCandidatesByDMA(graph.getNodes(), DMA);
	}

	/**
	 * Selects and returns a random node from a given list of nodes.
	 *
	 * @param nodes List of nodes from which to select randomly.
	 * @return A randomly selected node from the list, or {@code null} if the list
	 *         is empty (previously this threw on empty lists, forcing callers to
	 *         wrap every lookup in try/catch).
	 */
	public static NodeGraph selectRandomNode(List<NodeGraph> nodes) {
		return selectRandomNode(nodes, fallbackGenerator());
	}

	/**
	 * As {@link #selectRandomNode(List)}, drawing from the supplied generator.
	 *
	 * @param nodes  List of nodes from which to select randomly.
	 * @param random The generator to draw from.
	 * @return A randomly selected node from the list, or {@code null} if the list
	 *         is empty.
	 */
	public static NodeGraph selectRandomNode(List<NodeGraph> nodes, MersenneTwisterFast random) {
		if (nodes == null || nodes.isEmpty()) {
			return null;
		}
		return nodes.get(random.nextInt(nodes.size()));
	}
}
