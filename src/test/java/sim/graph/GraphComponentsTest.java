package sim.graph;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertNotSame;
import static org.junit.jupiter.api.Assertions.assertNull;
import static org.junit.jupiter.api.Assertions.assertSame;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.util.List;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.DisplayName;
import org.junit.jupiter.api.Test;
import sim.testing.Fixtures;

/** Behaviour of the two graph components, {@link NodeGraph} and {@link EdgeGraph}. */
class GraphComponentsTest {

  private Graph graph;

  @BeforeEach
  void buildGrid() {
    graph = Fixtures.grid(3, 3, 100.0);
  }

  @Test
  @DisplayName("a corner junction has two neighbours, an inner one four")
  void adjacencyFollowsTheGridDegree() {
    assertEquals(2, Fixtures.nodeAt(graph, 0.0, 0.0).getAdjacentNodes().size());
    assertEquals(4, Fixtures.nodeAt(graph, 100.0, 100.0).getAdjacentNodes().size());
    assertEquals(3, Fixtures.nodeAt(graph, 100.0, 0.0).getEdges().size());
  }

  @Test
  @DisplayName("getAdjacentNodes() hands back a copy")
  void adjacentNodesAreDefensivelyCopied() {
    NodeGraph node = Fixtures.nodeAt(graph, 100.0, 100.0);
    List<NodeGraph> first = node.getAdjacentNodes();
    first.clear();
    assertEquals(4, node.getAdjacentNodes().size());
    assertNotSame(first, node.getAdjacentNodes());
  }

  @Test
  void outDirectedEdgesMatchTheNodeDegree() {
    NodeGraph node = Fixtures.nodeAt(graph, 100.0, 100.0);
    assertEquals(4, node.getOutDirectedEdges().size());
  }

  @Test
  @DisplayName("getAdjacentRegion() is null unless the node is a gateway")
  void adjacentRegionOnlyAnswersForGateways() {
    NodeGraph node = Fixtures.nodeAt(graph, 100.0, 100.0);
    assertNull(node.getAdjacentRegion());

    node.setRegionID(1);
    for (NodeGraph neighbour : node.getAdjacentNodes()) {
      neighbour.setRegionID(2);
    }
    node.gateway = true;
    assertEquals(4, node.getAdjacentRegion().size());
  }

  @Test
  void nodeCarriesItsIdAndCentrality() {
    NodeGraph node = Fixtures.nodeAt(graph, 0.0, 0.0);
    node.setCentrality(12.5);
    assertEquals(12.5, node.getCentrality(), 0.0);
    assertTrue(graph.getNodeIDs().contains(node.getID()));
  }

  @Test
  @DisplayName("getOtherNode() walks across an edge, and is null off it")
  void otherNodeCrossesTheEdge() {
    NodeGraph origin = Fixtures.nodeAt(graph, 0.0, 0.0);
    NodeGraph neighbour = Fixtures.nodeAt(graph, 100.0, 0.0);
    NodeGraph stranger = Fixtures.nodeAt(graph, 200.0, 200.0);
    EdgeGraph edge = graph.getEdgeBetween(origin, neighbour);

    assertSame(neighbour, edge.getOtherNode(origin));
    assertSame(origin, edge.getOtherNode(neighbour));
    assertNull(edge.getOtherNode(stranger));
  }

  @Test
  void commonNodeFindsTheSharedJunction() {
    NodeGraph shared = Fixtures.nodeAt(graph, 100.0, 0.0);
    EdgeGraph west = graph.getEdgeBetween(Fixtures.nodeAt(graph, 0.0, 0.0), shared);
    EdgeGraph east = graph.getEdgeBetween(shared, Fixtures.nodeAt(graph, 200.0, 0.0));
    EdgeGraph elsewhere = graph.getEdgeBetween(Fixtures.nodeAt(graph, 0.0, 100.0),
        Fixtures.nodeAt(graph, 0.0, 200.0));

    assertSame(shared, west.getCommonNode(east));
    assertNull(west.getCommonNode(elsewhere));
  }

  @Test
  void edgeReportsItsLengthAndCentroid() {
    EdgeGraph edge = graph.getEdgeBetween(Fixtures.nodeAt(graph, 0.0, 0.0),
        Fixtures.nodeAt(graph, 100.0, 0.0));

    assertEquals(100.0, edge.getLength(), 1e-9);
    assertEquals(50.0, edge.getCoordsCentroid().x, 1e-9);
    assertEquals(0.0, edge.getCoordsCentroid().y, 1e-9);
    assertEquals(2, edge.getNodes().size());
  }

  @Test
  @DisplayName("getLine() hands back a copy of the geometry")
  void edgeLineIsDefensivelyCopied() {
    EdgeGraph edge = graph.getEdges().get(0);
    assertNotSame(edge.getLine(), edge.getLine());
    assertTrue(edge.getLine().equalsExact(edge.getLine()));
  }

  @Test
  void agentCountIsIncrementedDecrementedAndReset() {
    EdgeGraph edge = graph.getEdges().get(0);
    assertEquals(0, edge.getAgentCount());

    edge.incrementAgentCount();
    edge.incrementAgentCount();
    assertEquals(2, edge.getAgentCount());

    edge.decrementAgentCount();
    assertEquals(1, edge.getAgentCount());

    edge.resetAgentCount();
    assertEquals(0, edge.getAgentCount());
  }
}
