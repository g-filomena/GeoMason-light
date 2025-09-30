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
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;
import org.locationtech.jts.geom.Coordinate;
import org.locationtech.jts.geom.CoordinateArrays;
import org.locationtech.jts.geom.LineString;
import org.locationtech.jts.planargraph.DirectedEdge;
import sim.util.geo.Utilities;

/**
 * This class represents subgraphs derived from a graph. It establishes links between the components
 * of the parent graph and the corresponding child components, allowing for faster operations and
 * the creation of "regional" or "district" graphs.
 *
 * Agent navigation within a SubGraph is straightforward and can be easily retraced to the parent
 * graph.
 */
public class SubGraph extends Graph {

  private final SubGraphNodesMap subGraphNodesMap = new SubGraphNodesMap();
  private final SubGraphEdgesMap subGraphEdgesMap = new SubGraphEdgesMap();

  /**
   * Constructs a subgraph from a parent graph by copying a specified list of edges.
   *
   * @param parentGraph The graph from which to copy edges and nodes;
   */
  public SubGraph(Graph parentGraph) {
    edgesGraph = parentGraph.getEdges();
    nodesGraph = parentGraph.getNodes();
    for (NodeGraph node : nodesGraph) {
      node.setNeighbouringComponents();
    }
    generateAdjacencyMatrix();
  }

  /**
   * Constructs a subgraph from a parent graph by copying a specified list of edges.
   *
   * @param edges The list of edges to be included in the subgraph.
   */
  public SubGraph(List<EdgeGraph> edges) {
    for (EdgeGraph edge : edges) {
      addFromParentGraph(edge);
    }
    for (NodeGraph node : nodesGraph) {
      node.setNeighbouringComponents();
    }
    generateAdjacencyMatrix();
  }

  /**
   * Constructs an empty subgraph.
   */
  public SubGraph() {}

  /**
   * Adds an edge and its corresponding nodes to the subgraph from the parent graph.
   *
   * @param parentEdge The edge to be added to the subgraph, which is originally from the parent
   *        graph.
   */
  public void addFromParentGraph(EdgeGraph parentEdge) {

    final NodeGraph fromNode = parentEdge.getFromNode();
    final NodeGraph toNode = parentEdge.getToNode();
    final Coordinate fromNodeCoord = fromNode.getCoordinate();
    final Coordinate toNodeCoord = toNode.getCoordinate();

    final NodeGraph childFromNode = getNode(fromNodeCoord);
    final NodeGraph childToNode = getNode(toNodeCoord);

    final LineString line = (LineString) parentEdge.getMasonGeometry().geometry;
    final EdgeGraph childEdge = new EdgeGraph(line);
    final Coordinate[] coords = CoordinateArrays.removeRepeatedPoints(parentEdge.getCoordinates());

    final DirectedEdge de0 = new DirectedEdge(childFromNode, childToNode, coords[1], true);
    final DirectedEdge de1 =
        new DirectedEdge(childToNode, childFromNode, coords[coords.length - 2], false);
    childEdge.setDirectedEdges(de0, de1);

    childEdge.setNodes(childFromNode, childToNode);
    setAttributesChildEdge(childEdge, parentEdge);
    setAttributesChildNode(childFromNode, fromNode);
    setAttributesChildNode(childToNode, toNode);

    subGraphNodesMap.addChildParentPair(childFromNode, fromNode);
    subGraphNodesMap.addChildParentPair(childToNode, toNode);
    subGraphEdgesMap.add(childEdge, parentEdge);
    addEdge(childEdge);
  }

  /**
   * A mapping utility class for managing relationships between nodes in a subgraph.
   */
  private class SubGraphNodesMap {
    public Map<NodeGraph, NodeGraph> childParentMap = new HashMap<>();

    /**
     * Adds a mapping between a child node and its parent node in the subgraph.
     *
     * @param childNode The child node to be added to the subgraph.
     * @param parentNode The parent node corresponding to the child node.
     */
    private void addChildParentPair(NodeGraph childNode, NodeGraph parentNode) {
      childParentMap.put(childNode, parentNode);
    }

    /**
     * Retrieves the parent node corresponding to a given child node in the subgraph.
     *
     * @param childNode The child node for which the parent node is to be retrieved.
     * @return The parent node corresponding to the provided child node, or null if not found.
     */
    private NodeGraph findParent(NodeGraph childNode) {
      return childParentMap.get(childNode);
    }

    /**
     * Retrieves the child node corresponding to a given parent node in the subgraph.
     *
     * @param parentNode The parent node for which the child node is to be retrieved.
     * @return The child node corresponding to the provided parent node, or null if not found.
     */
    private NodeGraph findChild(NodeGraph parentNode) {
      return Utilities.getKeyFromValue(childParentMap, parentNode);
    }
  }

  /**
   * A mapping utility class for managing relationships between edges in a subgraph.
   */
  private class SubGraphEdgesMap {

    private final Map<EdgeGraph, EdgeGraph> childParentMap = new HashMap<>();

    /**
     * Adds a mapping between a child edge and its parent edge in the subgraph.
     *
     * @param edge The child edge to be added to the subgraph.
     * @param parentEdge The parent edge corresponding to the child edge.
     */
    public void add(EdgeGraph childEdge, EdgeGraph parentEdge) {
      childParentMap.put(childEdge, parentEdge);
    }

    /**
     * Retrieves the parent edge corresponding to a given child edge in the subgraph.
     *
     * @param edgeSubGraph The child edge for which the parent edge is to be retrieved.
     * @return The parent edge corresponding to the provided child edge, or null if not found.
     */
    private EdgeGraph findParent(EdgeGraph childEdge) {
      return childParentMap.get(childEdge);
    }

    /**
     * Retrieves the child edge corresponding to a given parent edge in the subgraph.
     *
     * @param edgeSubGraph The parent edge for which the child edge is to be retrieved.
     * @return The child edge corresponding to the provided parent edge, or null if not found.
     */
    private EdgeGraph findChild(EdgeGraph parentGraph) {
      return Utilities.getKeyFromValue(childParentMap, parentGraph);
    }
  }

  public void updateSubGraph() {
    for (NodeGraph node : nodesGraph) {
      node.setNeighbouringComponents();
    }
    generateAdjacencyMatrix();
  }

  /**
   * Retrieves the parent node corresponding to the provided child node.
   *
   * @param childNode The child node for which the parent node is to be retrieved.
   * @return The parent node corresponding to the provided child node, or null if not found.
   */
  public NodeGraph getParentNode(NodeGraph childNode) {
    return this.subGraphNodesMap.findParent(childNode);
  }

  /**
   * Retrieves a List of all parent nodes within the current subgraph.
   *
   * @return A List of all parent nodes present in the current subgraph.
   */
  public List<NodeGraph> getParentNodes() {
    final List<NodeGraph> parentNodes = new ArrayList<>();
    parentNodes.addAll(subGraphNodesMap.childParentMap.values());
    return parentNodes;
  }

  /**
   * Retrieves the child node corresponding to the provided parent node.
   *
   * @param parentNode The parent node for which the child node is to be retrieved.
   * @return The child node corresponding to the provided parent node.
   */
  public NodeGraph getChildNode(NodeGraph parentNode) {

    final NodeGraph childNode = subGraphNodesMap.findChild(parentNode);
    return childNode;
  }

  /**
   * Retrieves the parent edge corresponding to the provided child edge.
   *
   * @param childEdge The child edge for which the parent edge is to be retrieved.
   * @return The parent edge corresponding to the provided child edge, or null if not found.
   */
  public EdgeGraph getParentEdge(EdgeGraph childEdge) {
    return subGraphEdgesMap.findParent(childEdge);
  }

  /**
   * Retrieves a List of all parent edges within the current subgraph.
   *
   * @return A List of all parent edges present in the current subgraph.
   */
  public List<EdgeGraph> getParentEdges() {
    final List<EdgeGraph> parentEdges = new ArrayList<>();
    parentEdges.addAll(subGraphEdgesMap.childParentMap.values());
    return parentEdges;
  }

  /**
   * Retrieves the child edge corresponding to the provided parent edge.
   *
   * @param parentEdge The parent edge for which the child edge is to be retrieved.
   * @return The child edge corresponding to the provided parent edge.
   */
  public EdgeGraph getChildEdge(EdgeGraph parentEdge) {

    final EdgeGraph childEdge = subGraphEdgesMap.findChild(parentEdge);
    return childEdge;
  }

  /**
   * Retrieves a List of parent nodes corresponding to the provided list of child nodes.
   *
   * @param childNodes The List of child nodes for which parent nodes are to be retrieved.
   * @return A List of parent nodes corresponding to the provided child nodes.
   */
  public List<NodeGraph> getParentNodes(List<NodeGraph> childNodes) {
    final List<NodeGraph> parentNodes = new ArrayList<>();
    for (final NodeGraph child : childNodes) {
      final NodeGraph parent = subGraphNodesMap.findParent(child);
      if (parent != null) {
        parentNodes.add(parent);
      }
    }
    return parentNodes;
  }

  /**
   * Retrieves a List of child nodes corresponding to the provided list of parent nodes.
   *
   * @param parentNodes The List of parent nodes for which child nodes are to be retrieved.
   * @return A List of child nodes corresponding to the provided parent nodes.
   */
  public List<NodeGraph> getChildNodes(List<NodeGraph> parentNodes) {
    final List<NodeGraph> childNodes = new ArrayList<>();
    for (final NodeGraph parent : parentNodes) {
      final NodeGraph child = subGraphNodesMap.findChild(parent);
      if (child != null) {
        childNodes.add(child);
      }
    }
    return childNodes;
  }

  /**
   * Retrieves the parent edges corresponding to the provided child edges.
   *
   * @param childEdges The List of child edges for which parent edges are to be retrieved.
   * @return A List of parent edges corresponding to the provided child edges.
   */
  public List<EdgeGraph> getParentEdges(List<EdgeGraph> childEdges) {
    final List<EdgeGraph> parentEdges = new ArrayList<>();
    for (final EdgeGraph child : childEdges) {
      final EdgeGraph parent = subGraphEdgesMap.findParent(child);
      if (parent != null) {
        parentEdges.add(parent);
      }
    }
    return parentEdges;
  }

  /**
   * Retrieves the child edges corresponding to the provided parent edges.
   *
   * @param parentEdges The List of parent edges for which child edges are to be retrieved.
   * @return A List of child edges corresponding to the provided parent edges.
   */
  public List<EdgeGraph> getChildEdges(List<EdgeGraph> parentEdges) {
    List<EdgeGraph> childEdges = new ArrayList<>();
    for (EdgeGraph parent : parentEdges) {
      EdgeGraph child = subGraphEdgesMap.findChild(parent);
      if (child != null) {
        childEdges.add(child);
      }
    }
    return childEdges;
  }

  /**
   * Sets attributes of a child node based on the attributes of the parentNode.
   *
   * This method is used to propagate certain attributes from a parent node to its child node.
   *
   * @param childNode The child node to which attributes will be set.
   * @param parentNode The parent node from which attributes will be copied.
   */
  public void setAttributesChildNode(NodeGraph childNode, NodeGraph parentNode) {

    childNode.nodeID = parentNode.getID();
    childNode.attributes = parentNode.attributes;
    // since this contains further attributes
    childNode.masonGeometry = parentNode.masonGeometry;
    childNode.regionID = parentNode.regionID;

    // for primalGraph only
    // childNode.visibleBuildings2d = parentNode.visibleBuildings2d;
    childNode.adjacentBuildings = parentNode.adjacentBuildings;
    childNode.visibleBuildings3d = parentNode.visibleBuildings3d;
    // TODO
    childNode.adjacentRegions = parentNode.adjacentRegions;
    childNode.DMA = parentNode.DMA;

    // for dualGraph only
    childNode.primalEdge = parentNode.primalEdge;
  }

  /**
   * Sets attributes of a child edge based on the attributes of a parent edge.
   *
   * This method is used to propagate certain attributes from a parent edge to its child edge, such
   * as the edge ID and attributes related to dual edges.
   *
   * @param childEdge The child edge to which attributes will be set.
   * @param parentEdge The parent edge from which attributes will be copied.
   */
  public void setAttributesChildEdge(EdgeGraph childEdge, EdgeGraph parentEdge) {

    childEdge.edgeID = parentEdge.getID();
    childEdge.attributes = parentEdge.attributes;
    childEdge.regionID = parentEdge.regionID;
    childEdge.masonGeometry = parentEdge.masonGeometry;

    // for primal edges:
    childEdge.dualNode = parentEdge.getDualNode();

    // for dual edges:
    childEdge.deflectionAngle = parentEdge.deflectionAngle;
  }

  /**
   * Generates the centrality map for the nodes of the current subgraph nodes from the centrality
   * values of their corresponding parent nodes in the parent graph. The resulting centrality values
   * for nodes in the subgraph are sorted in descending order and stored as the subgraph's
   * centrality map.
   */
  public void generateSubGraphCentralityMap() {

    Map<NodeGraph, Double> centralityMap = new LinkedHashMap<>();
    List<NodeGraph> nodes = nodesGraph;
    for (NodeGraph node : nodes) {
      NodeGraph parentNode = getParentNode(node);
      if (parentNode.centrality == Double.MAX_VALUE) {
        return;
      }
      centralityMap.put(node, parentNode.centrality);
    }
    this.centralityMap = Utilities.sortByValue(centralityMap, false);
  }

  /**
   * Calculates and returns salient nodes within a subgraph based on a specified percentile of
   * centrality values. This method identifies nodes with centrality values equal to or higher than
   * the specified percentile within the subgraph, then maps those nodes to their parent nodes in
   * the parent graph and returns the resulting subgraph's salient nodes.
   *
   * @param percentile The percentile value used to determine salient nodes.
   * @return A mapping of salient nodes within the subgraph to their centrality values, or null if
   *         no salient nodes are found in the subgraph.
   */
  public Map<NodeGraph, Double> getSubGraphSalientNodes(double percentile) {

    int position;
    position = (int) (centralityMap.size() * percentile);
    final double boundary = new ArrayList<>(centralityMap.values()).get(position);

    final Map<NodeGraph, Double> filteredMap =
        centralityMap.entrySet().stream().filter(entry -> entry.getValue() >= boundary)
            .collect(Collectors.toMap(entry -> entry.getKey(), entry -> entry.getValue()));

    if (filteredMap.isEmpty() || filteredMap == null) {
      return null;
    }
    Map<NodeGraph, Double> parentMap = new HashMap<>();

    for (Map.Entry<NodeGraph, Double> entry : filteredMap.entrySet()) {
      NodeGraph parentNode = this.getParentNode(entry.getKey());
      parentMap.put(parentNode, entry.getValue());
    }
    return parentMap;
  }

}
