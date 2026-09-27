/*
 * Copyright (c) 2023 Gabriele Filomena University of Liverpool, UK
 *
 * This program is free software: it can redistributed and/or modified under the terms of the GNU
 * General Public License 3.0 as published by the Free Software Foundation.
 *
 * See the file "LICENSE" for more information
 */
package sim.graph;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertNotNull;

import java.util.Arrays;
import java.util.Collections;
import java.util.Map;
import org.junit.jupiter.api.DisplayName;
import org.junit.jupiter.api.Test;
import sim.testing.Fixtures;

class NodeGraphTest {

  /** Gives every edge a dual node at its centroid, all in region 1. */
  private static void attachDualNodes(Graph graph) {
    for (EdgeGraph edge : graph.getEdges()) {
      NodeGraph dualNode = new NodeGraph(edge.getCoordsCentroid());
      dualNode.setPrimalEdge(edge);
      edge.setDualNode(dualNode);
      edge.setRegionID(1);
    }
  }

  @Test
  @DisplayName("region-based dual nodes can be asked of a node other than the origin")
  void dualNodesOfAnIntermediateNode() {
    Graph graph = Fixtures.grid(3, 3, 100.0);
    attachDualNodes(graph);
    NodeGraph origin = Fixtures.nodeAt(graph, 0, 0);
    NodeGraph middle = Fixtures.nodeAt(graph, 100, 100);
    NodeGraph destination = Fixtures.nodeAt(graph, 200, 200);

    Map<NodeGraph, Double> dualNodes = middle.getDualNodes(origin, destination, true, null);

    assertEquals(4, dualNodes.size());
    assertNotNull(middle.getDualNode(origin, destination, true, null));
  }

  @Test
  @DisplayName("an edge without a dual node is not offered as a departure")
  void edgesWithoutADualNodeAreSkipped() {
    Graph graph = Fixtures.path(3, 100.0);
    EdgeGraph first = graph.getEdgeBetween(Fixtures.nodeAt(graph, 0, 0),
        Fixtures.nodeAt(graph, 100, 0));
    NodeGraph dualNode = new NodeGraph(first.getCoordsCentroid());
    dualNode.setPrimalEdge(first);
    first.setDualNode(dualNode);

    Map<NodeGraph, Double> dualNodes = Fixtures.nodeAt(graph, 100, 0)
        .getDualNodes(Fixtures.nodeAt(graph, 0, 0), Fixtures.nodeAt(graph, 200, 0), false, null);

    assertEquals(Collections.singleton(dualNode), dualNodes.keySet());
    assertFalse(dualNodes.containsKey(null));
  }

  @Test
  @DisplayName("asking for a gateway's adjacent regions twice does not duplicate its entries")
  void adjacentRegionEntriesAreRebuilt() {
    Graph graph = Fixtures.path(2, 100.0);
    NodeGraph gateway = Fixtures.nodeAt(graph, 0, 0);
    NodeGraph other = Fixtures.nodeAt(graph, 100, 0);
    gateway.gateway = true;
    other.setRegionID(3);

    gateway.getAdjacentRegion();

    assertEquals(Arrays.asList(3), gateway.getAdjacentRegion());
    assertEquals(Arrays.asList(other), gateway.adjacentRegionEntries);
  }

  @Test
  @DisplayName("setting an edge's nodes again replaces them")
  void setNodesReplaces() {
    Graph graph = Fixtures.path(2, 100.0);
    EdgeGraph edge = graph.getEdges().get(0);

    edge.setNodes(edge.getFromNode(), edge.getToNode());

    assertEquals(2, edge.getNodes().size());
  }
}
