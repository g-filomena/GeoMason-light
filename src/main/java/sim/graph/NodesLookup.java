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
import java.util.Random;
import java.util.stream.Collectors;
import sim.field.geo.VectorLayer;
import sim.util.geo.MasonGeometry;

/**
 * A class containing functions to identify random nodes within a graph, based on specific criteria
 * such as location, region, and distance. These methods are useful for generating origin and
 * destination nodes for agent-based simulations in spatial models.
 */
public class NodesLookup {

  private static final Random random = new Random();
  final static double PERCENTILE_DECREASE = 0.05;
  final static double EXPANSION_FACTOR = 1.10;
  final static double RADIUS_THRESHOLD = 2.00;

  final static double MIN_DISTANCE = 100;
  final static double INITIAL_TOLERANCE = 50;
  final static double TOLERANCE_INCREMENT = 50;

  final static double DISTANCE_MULTIPLIER = 1.50;

  /**
   * Returns a randomly selected node from a graph. This method is useful for generating random
   * starting or ending points in simulations.
   *
   * @param graph The graph from which a random node is to be selected.
   * @return A randomly chosen NodeGraph object from the graph.
   */
  public static NodeGraph randomNode(Graph graph) {
    final List<NodeGraph> candidates = new ArrayList<>(graph.nodesGraph);
    return selectRandomNode(candidates);
  }

  /**
   * Returns a randomly selected node from a specified list of NodedGraph elements.
   *
   * @param nodes The list of nodes.
   * @return A NodeGraph object randomly selected from the given list.
   */
  public static NodeGraph randomNodeFromList(List<NodeGraph> nodes) {

    int randomInt = random.nextInt(nodes.size());
    return nodes.get(randomInt);
  }

  /**
   * Returns a randomly selected node from a specified list of node geometries within a graph. This
   * method allows for the random selection of nodes that meet certain spatial criteria represented
   * by the node geometries.
   *
   * @param graph The graph containing the nodes.
   * @param nodesGeometries A List of MasonGeometry objects representing specific node locations.
   * @return A NodeGraph object randomly selected from the specified geometries.
   */
  public static NodeGraph randomNodeFromGeometriest(Graph graph,
      List<MasonGeometry> nodesGeometries) {
    Integer randomInt = random.nextInt(nodesGeometries.size());
    MasonGeometry geoNode = nodesGeometries.get(randomInt);
    return graph.findNode(geoNode.geometry.getCoordinate());
  }

  /**
   * Returns a randomly selected node from a graph within a specified radius from an origin node and
   * outside the origin node's region. This method is used to generate nodes that are geographically
   * spread yet regionally distinct from the origin.
   *
   * @param graph The graph from which to select the node.
   * @param originNode The origin node serving as the center of the search radius.
   * @param radius The radius within which to search for a suitable node.
   * @return A randomly selected NodeGraph object that meets the specified criteria.
   */
  public static NodeGraph randomNodeRegion(Graph graph, NodeGraph originNode, double radius) {

    final MasonGeometry originNodeGeometry = originNode.masonGeometry;
    double expandingRadius = radius;

    while (true) {
      if (expandingRadius >= radius * RADIUS_THRESHOLD) {
        return null;
      }

      List<MasonGeometry> spatialFilter =
          graph.junctions.featuresWithinDistance(originNodeGeometry.geometry, expandingRadius);
      if (spatialFilter.isEmpty()) {
        expandingRadius *= EXPANSION_FACTOR;
        continue;
      }

      List<MasonGeometry> regionFilter = spatialFilter.stream()
          .filter(
              geo -> Integer.parseInt(geo.getStringAttribute("district")) == originNode.regionID)
          .collect(Collectors.toList());

      if (regionFilter.isEmpty()) {
        expandingRadius *= EXPANSION_FACTOR;
        continue;
      }

      MasonGeometry nodeGeometry = regionFilter.get(random.nextInt(regionFilter.size()));
      NodeGraph node = graph.findNode(nodeGeometry.geometry.getCoordinate());

      if (node != null) {
        return node;
      }

      expandingRadius *= EXPANSION_FACTOR;
    }
  }

  /**
   * Returns a randomly selected node from a graph whose distance from the origin node matches one
   * of the provided distances.
   *
   * @param graph The input graph.
   * @param junctions The vector layer representing junctions.
   * @param originNode The origin node.
   * @param distances The list of possible distances used to identify the node.
   * @return A randomly selected node that matches the specified distance criteria.
   */
  public static NodeGraph randomNodeFromDistancesSet(Graph graph, VectorLayer junctions,
      NodeGraph originNode, List<Float> distances) {

    // Select a random distance from the list
    double distance = distances.get(random.nextInt(distances.size()));
    if (distance < MIN_DISTANCE) {
      distance = MIN_DISTANCE;
    }

    NodeGraph node = null;
    List<NodeGraph> candidates = new ArrayList<>();
    double tolerance = INITIAL_TOLERANCE;

    while (candidates.isEmpty()) {
      double lowerLimit = distance - tolerance;
      double upperLimit = distance + tolerance;
      candidates = getNodesBetweenDistanceInterval(graph, originNode, lowerLimit, upperLimit);
      tolerance += TOLERANCE_INCREMENT;
    }

    while (node == null || node.getID() == originNode.getID()) {
      node = candidates.get(random.nextInt(candidates.size()));
      if (graph.getEdgeBetween(originNode, node) != null) {
        node = null;
      }
    }

    return node;
  }

  /**
   * Returns a List of nodes in the graph whose distance from the given node falls within the
   * specified range.
   *
   * @param graph The input graph.
   * @param node The reference node.
   * @param lowerLimit The minimum distance from the reference node.
   * @param upperLimit The maximum distance from the reference node.
   * @return A List of nodes that meet the distance criteria.
   */
  public static List<NodeGraph> getNodesBetweenDistanceInterval(Graph graph, NodeGraph node,
      double lowerLimit, double upperLimit) {
    MasonGeometry originGeometry = node.masonGeometry;
    List<MasonGeometry> containedGeometries =
        graph.junctions.featuresBetweenLimits(originGeometry.geometry, lowerLimit, upperLimit);
    return containedGeometries.stream()
        .map(masonGeometry -> graph.findNode(masonGeometry.geometry.getCoordinate()))
        .filter(potentialNode -> !node.equals(potentialNode)).collect(Collectors.toList());
  }

  /**
   * Returns a randomly selected node from a graph whose distance from an origin node falls within a
   * specified range. This method allows for the selection of nodes based on proximity constraints.
   *
   * @param graph The graph to search within.
   * @param originNode The origin node to measure distances from.
   * @param lowerLimit The minimum distance from the origin node.
   * @param upperLimit The maximum distance from the origin node.
   * @return A NodeGraph object randomly selected within the specified distance range.
   */
  public static NodeGraph randomNodeBetweenDistanceInterval(Graph graph, NodeGraph originNode,
      double lowerLimit, double upperLimit) {

    final List<NodeGraph> candidates =
        getNodesBetweenDistanceInterval(graph, originNode, lowerLimit, upperLimit);
    return selectRandomNode(candidates);
  }

  /**
   * Returns a List of nodes whose distance from the given node falls within the specified range and
   * belong to a different region.
   *
   * @param graph The input graph.
   * @param node The reference node.
   * @param lowerLimit The minimum distance from the reference node.
   * @param upperLimit The maximum distance from the reference node.
   * @return A List of nodes that meet the distance and region criteria.
   */
  public static List<NodeGraph> getNodesBetweenDistanceIntervalRegion(Graph graph, NodeGraph node,
      double lowerLimit, double upperLimit) {
    List<NodeGraph> containedNodes = new ArrayList<>();
    containedNodes = getNodesBetweenDistanceInterval(graph, node, lowerLimit, upperLimit);
    return graph.nodesInRegion(containedNodes, node.regionID);
  }

  /**
   * Returns a randomly selected node from a graph whose distance from an origin node falls within a
   * specified range and belongs to a different region than the origin. This method is useful for
   * simulating regional travel or migration patterns.
   *
   * @param graph The graph to search within.
   * @param originNode The origin node to measure distances from.
   * @param lowerLimit The minimum distance from the origin node.
   * @param upperLimit The maximum distance from the origin node.
   * @return A NodeGraph object randomly selected within the specified distance range and different
   *         region.
   */
  public static NodeGraph randomNodeBetweenDistanceIntervalRegion(Graph graph, NodeGraph originNode,
      double lowerLimit, double upperLimit) {

    List<NodeGraph> candidates =
        getNodesBetweenDistanceIntervalRegion(graph, originNode, lowerLimit, upperLimit);
    return selectRandomNode(candidates);
  }

  /**
   * Returns a randomly selected node whose distance from the origin node falls within the specified
   * range and has centrality values above or equal to a specified percentile.
   *
   * @param graph The input graph.
   * @param originNode The origin node.
   * @param lowerLimit The minimum distance from the origin node.
   * @param upperLimit The maximum distance from the origin node.
   * @param percentile The percentile used as a threshold for centrality values.
   * @return A randomly selected node that satisfies the specified distance and centrality criteria.
   */
  public static NodeGraph randomSalientNodeBetweenDistanceInterval(Graph graph,
      NodeGraph originNode, double lowerLimit, double upperLimit, double percentile) {

    NodeGraph node = null;

    while (node == null) {
      List<NodeGraph> candidates = graph.getSalientNodesBetweenDistanceInterval(originNode,
          lowerLimit, upperLimit, percentile);

      if (candidates.isEmpty()) {
        return null; // Return null if no candidates are found
      }

      node = candidates.get(random.nextInt(candidates.size()));

      // Reduce percentile for the next iteration if node is not found
      percentile -= PERCENTILE_DECREASE;
      if (percentile <= 0.0) {
        return null; // Return null if percentile reaches 0 or below
      }
    }
    return node;
  }

  /**
   * Returns a random node from the graph whose distance from the origin node falls within a
   * specified range and belongs to a specific category based on the DMA label. This method is
   * useful for selecting nodes that are relevant to specific activities or functions like living,
   * working, or visiting.
   *
   * DMA = Urban Density (D), Mix(M) and Access (A).
   *
   * @param graph The graph to search within.
   * @param originNode The origin node to measure distances from.
   * @param lowerLimit The minimum distance from the origin node.
   * @param upperLimit The maximum distance from the origin node.
   * @param DMA The category label (e.g., "live," "work," "visit") or "random" for any category.
   * @return A NodeGraph object randomly selected based on the specified distance and category
   *         criteria.
   */
  public static NodeGraph randomNodeBetweenDistanceIntervalDMA(Graph graph, NodeGraph originNode,
      double lowerLimit, double upperLimit, String DMA) {

    NodeGraph node = null;

    while (node == null) {
      List<NodeGraph> candidates =
          getNodesBetweenDistanceInterval(graph, originNode, lowerLimit, upperLimit);
      List<NodeGraph> candidatesDMA = getCandidatesByDMA(candidates, DMA);

      if (candidatesDMA.isEmpty()) {
        if (upperLimit > upperLimit * DISTANCE_MULTIPLIER) {
          node = candidates.get(random.nextInt(candidates.size()));
          break;
        }
        upperLimit += INITIAL_TOLERANCE;
      } else {
        node = candidatesDMA.get(random.nextInt(candidatesDMA.size()));
      }
    }
    return node;
  }

  /**
   * Returns a randomly selected node from the graph that belongs to a specified category ("live,"
   * "work," "visit") based on the provided DMA label.
   *
   * DMA = Urban Density (D), Mix(M) and Access (A).
   *
   * @param graph The input graph.
   * @param DMA The desired node category ("live," "work," "visit") or "random".
   * @return A randomly selected node from the specified category based on the DMA.
   */
  public static NodeGraph randomNodeDMA(Graph graph, String DMA) {
    List<NodeGraph> candidates = graph.getNodes();
    List<NodeGraph> candidatesDMA = getCandidatesByDMA(candidates, DMA);
    return selectRandomNode(candidatesDMA);
  }

  /**
   * Filters a list of nodes based on the specified DMA. - If DMA is "random", returns nodes whose
   * DMA is either "work", "visit", or "live". - If DMA is "workOrVisit", returns nodes whose DMA is
   * either "work" or "visit". - Otherwise, returns nodes whose DMA exactly matches the given DMA
   * string.
   *
   * DMA = Urban Density (D), Mix(M) and Access (A).
   *
   * @param nodes List of nodes to filter.
   * @param DMA The DMA criterion used for filtering ("random", "workOrVisit", or specific DMA).
   * @return A filtered list of nodes matching the DMA criteria.
   */
  public static List<NodeGraph> getCandidatesByDMA(List<NodeGraph> nodes, String DMA) {
    if (DMA.equals("random")) {
      return nodes.stream().filter(
          node -> node.DMA.equals("work") || node.DMA.equals("visit") || node.DMA.equals("live"))
          .collect(Collectors.toList());
    } else if (DMA.equals("workOrVisit")) {
      return nodes.stream().filter(node -> node.DMA.equals("work") || node.DMA.equals("visit"))
          .collect(Collectors.toList());
    } else {
      return nodes.stream().filter(node -> node.DMA.equals(DMA)).collect(Collectors.toList());
    }
  }

  /**
   * Retrieves nodes from a graph based on the specified DMA.
   *
   * DMA = Urban Density (D), Mix(M) and Access (A).
   *
   * @param graph The graph containing the nodes to filter.
   * @param DMA The DMA criterion used for filtering ("random", "workOrVisit", or specific DMA).
   * @return A filtered list of nodes matching the DMA criteria.
   */
  public static List<NodeGraph> getNodesByDMA(Graph graph, String DMA) {
    return getCandidatesByDMA(graph.getNodes(), DMA);
  }

  /**
   * Selects and returns a random node from a given list of nodes.
   *
   * @param nodes List of nodes from which to select randomly.
   * @return A randomly selected node from the list.
   */
  public static NodeGraph selectRandomNode(List<NodeGraph> nodes) {
    return nodes.get(random.nextInt(nodes.size()));
  }
}
