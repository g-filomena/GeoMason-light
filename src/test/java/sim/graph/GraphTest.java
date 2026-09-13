package sim.graph;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertNotNull;
import static org.junit.jupiter.api.Assertions.assertNull;
import static org.junit.jupiter.api.Assertions.assertSame;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.DisplayName;
import org.junit.jupiter.api.Test;
import org.locationtech.jts.geom.Coordinate;
import org.locationtech.jts.geom.Polygon;
import org.locationtech.jts.planargraph.DirectedEdge;
import sim.testing.Fixtures;

class GraphTest {

  private Graph graph;

  @BeforeEach
  void buildGrid() {
    graph = Fixtures.grid(3, 3, 100.0);
  }

  @Test
  @DisplayName("a 3x3 grid has 9 junctions and 12 segments")
  void gridHasTheExpectedNodesAndEdges() {
    assertEquals(9, graph.getNodes().size());
    assertEquals(12, graph.getEdges().size());
    assertEquals(9, graph.getNodeIDs().size());
    assertEquals(12, graph.getEdgeIDs().size());
  }

  @Test
  void findNodeLocatesAJunctionByCoordinate() {
    assertNotNull(graph.findNode(new Coordinate(100.0, 100.0)));
    assertNull(graph.findNode(new Coordinate(1000.0, 1000.0)));
  }

  @Test
  void findNodeAcceptsANodeAsALookupKey() {
    NodeGraph node = Fixtures.nodeAt(graph, 100.0, 0.0);
    assertSame(node, graph.findNode(node));
  }

  @Test
  @DisplayName("getEdgeBetween() is symmetric and null for unconnected junctions")
  void edgeLookupIsSymmetric() {
    NodeGraph origin = Fixtures.nodeAt(graph, 0.0, 0.0);
    NodeGraph neighbour = Fixtures.nodeAt(graph, 100.0, 0.0);
    NodeGraph diagonal = Fixtures.nodeAt(graph, 100.0, 100.0);

    EdgeGraph edge = graph.getEdgeBetween(origin, neighbour);
    assertNotNull(edge);
    assertSame(edge, graph.getEdgeBetween(neighbour, origin));
    assertNull(graph.getEdgeBetween(origin, diagonal));
  }

  @Test
  void directedEdgeLookupFollowsTheRequestedDirection() {
    NodeGraph origin = Fixtures.nodeAt(graph, 0.0, 0.0);
    NodeGraph neighbour = Fixtures.nodeAt(graph, 100.0, 0.0);

    DirectedEdge outgoing = graph.getDirectedEdgeBetween(origin, neighbour);
    DirectedEdge incoming = graph.getDirectedEdgeBetween(neighbour, origin);

    assertNotNull(outgoing);
    assertNotNull(incoming);
    assertSame(origin, outgoing.getFromNode());
    assertSame(neighbour, outgoing.getToNode());
    assertSame(outgoing.getSym(), incoming);
  }

  @Test
  @DisplayName("the indexed and the linear containment queries agree")
  void nodesWithinPolygonMatchesContainedNodes() {
    Polygon polygon = (Polygon) Fixtures.FACTORY.toGeometry(
        new org.locationtech.jts.geom.Envelope(-10.0, 110.0, -10.0, 110.0));

    List<NodeGraph> indexed = graph.getNodesWithinPolygon(polygon);
    List<NodeGraph> scanned = graph.getContainedNodes(polygon);

    assertEquals(4, indexed.size());
    assertEquals(scanned.size(), indexed.size());
    assertTrue(indexed.containsAll(scanned));
  }

  @Test
  void containedEdgesPicksUpSegmentsInsideAGeometry() {
    Polygon polygon = (Polygon) Fixtures.FACTORY.toGeometry(
        new org.locationtech.jts.geom.Envelope(-10.0, 110.0, -10.0, 110.0));
    // The four segments of the bottom-left block are inside; those leaving it are not.
    assertEquals(4, graph.getContainedEdges(polygon).size());
  }

  @Test
  void nodesInRegionDropsTheGivenRegion() {
    List<NodeGraph> nodes = graph.getNodes();
    nodes.get(0).setRegionID(1);
    nodes.get(1).setRegionID(1);
    for (int i = 2; i < nodes.size(); i++) {
      nodes.get(i).setRegionID(2);
    }
    assertEquals(nodes.size() - 2, graph.nodesInRegion(nodes, 1).size());
  }

  @Test
  @DisplayName("getSalientNodes() keeps the top slice of the centrality map")
  void salientNodesKeepTheTopSlice() {
    List<NodeGraph> nodes = graph.getNodes();
    for (int i = 0; i < nodes.size(); i++) {
      nodes.get(i).setCentrality(i);
    }
    graph.generateCentralityMap();

    Map<NodeGraph, Double> salient = graph.getSalientNodes(0.75);
    // 9 nodes, position 6 in the ascending order, so centralities 6, 7 and 8 survive.
    assertEquals(3, salient.size());
    for (Double centrality : salient.values()) {
      assertTrue(centrality >= 6.0);
    }
  }

  @Test
  @DisplayName("a percentile of 1 keeps the single most central node instead of throwing")
  void salientNodesHandleTheTopPercentile() {
    List<NodeGraph> nodes = graph.getNodes();
    for (int i = 0; i < nodes.size(); i++) {
      nodes.get(i).setCentrality(i);
    }
    graph.generateCentralityMap();

    Map<NodeGraph, Double> salient = graph.getSalientNodes(1.0);
    assertEquals(1, salient.size());
    assertTrue(salient.containsValue(8.0));
  }

  @Test
  void salientNodesOfAnEmptyCentralityMapAreEmpty() {
    assertTrue(new Graph().getSalientNodes(0.5).isEmpty());
  }

  @Test
  @DisplayName("the spatially filtered salient lookup also survives a percentile of 1")
  void salientNodesWithinSpaceHandleTheTopPercentile() {
    List<NodeGraph> nodes = graph.getNodes();
    for (int i = 0; i < nodes.size(); i++) {
      nodes.get(i).setCentrality(i);
    }
    graph.generateCentralityMap();

    Map<NodeGraph, Double> salient = graph.getSalientNodesWithinSpace(
        Fixtures.nodeAt(graph, 0.0, 0.0), Fixtures.nodeAt(graph, 200.0, 200.0), 1.0);
    assertEquals(1, salient.size());
  }

  @Test
  void filterCentralityMapRetainsOnlyTheListedNodes() {
    List<NodeGraph> nodes = graph.getNodes();
    Map<NodeGraph, Double> centrality = new LinkedHashMap<>();
    for (int i = 0; i < nodes.size(); i++) {
      centrality.put(nodes.get(i), (double) i);
    }

    List<NodeGraph> keep = new ArrayList<>(Arrays.asList(nodes.get(0), nodes.get(3)));
    Map<NodeGraph, Double> filtered = Graph.filterCentralityMap(centrality, keep);

    assertEquals(2, filtered.size());
    assertTrue(filtered.keySet().containsAll(keep));
  }

  @Test
  void edgesInNodesSpaceCoversTheSegmentsBetweenTwoJunctions() {
    NodeGraph origin = Fixtures.nodeAt(graph, 0.0, 0.0);
    NodeGraph destination = Fixtures.nodeAt(graph, 200.0, 200.0);

    List<EdgeGraph> edges = graph.edgesInNodesSpace(origin, destination);
    assertFalse(edges.isEmpty());
    assertTrue(graph.getEdges().containsAll(edges));
  }

  @Test
  @DisplayName("getNode() reuses the junction already at a coordinate")
  void getNodeIsIdempotentForAKnownCoordinate() {
    int before = graph.getNodes().size();
    NodeGraph existing = graph.getNode(new Coordinate(100.0, 100.0));
    assertSame(Fixtures.nodeAt(graph, 100.0, 100.0), existing);
    assertEquals(before, graph.getNodes().size());
  }
}
