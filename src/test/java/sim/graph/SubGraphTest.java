package sim.graph;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertNotNull;
import static org.junit.jupiter.api.Assertions.assertNotSame;
import static org.junit.jupiter.api.Assertions.assertNull;
import static org.junit.jupiter.api.Assertions.assertSame;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.DisplayName;
import org.junit.jupiter.api.Test;
import sim.testing.Fixtures;

class SubGraphTest {

  private Graph parent;
  private List<EdgeGraph> block;

  @BeforeEach
  void buildParentAndBlock() {
    parent = Fixtures.grid(3, 3, 100.0);

    // The four segments around the bottom-left block of the grid.
    NodeGraph southWest = Fixtures.nodeAt(parent, 0.0, 0.0);
    NodeGraph southEast = Fixtures.nodeAt(parent, 100.0, 0.0);
    NodeGraph northEast = Fixtures.nodeAt(parent, 100.0, 100.0);
    NodeGraph northWest = Fixtures.nodeAt(parent, 0.0, 100.0);

    block = new ArrayList<>();
    block.add(parent.getEdgeBetween(southWest, southEast));
    block.add(parent.getEdgeBetween(southEast, northEast));
    block.add(parent.getEdgeBetween(northEast, northWest));
    block.add(parent.getEdgeBetween(northWest, southWest));
  }

  @Test
  @DisplayName("a subgraph holds its own copies of the selected edges and their nodes")
  void subGraphCopiesTheSelectedEdges() {
    SubGraph sub = new SubGraph(block);

    assertEquals(4, sub.getEdges().size());
    assertEquals(4, sub.getNodes().size());
    for (EdgeGraph childEdge : sub.getEdges()) {
      assertTrue(block.contains(sub.getParentEdge(childEdge)));
    }
    for (NodeGraph childNode : sub.getNodes()) {
      assertNotSame(childNode, sub.getParentNode(childNode));
      assertNotNull(sub.getParentNode(childNode));
    }
  }

  @Test
  @DisplayName("the child and parent lookups are inverses of each other")
  void childAndParentLookupsRoundTrip() {
    SubGraph sub = new SubGraph(block);

    for (EdgeGraph parentEdge : block) {
      EdgeGraph childEdge = sub.getChildEdge(parentEdge);
      assertNotNull(childEdge);
      assertSame(parentEdge, sub.getParentEdge(childEdge));
      assertEquals(parentEdge.getID(), childEdge.getID());
      assertEquals(parentEdge.getLength(), childEdge.getLength(), 1e-9);
    }

    for (NodeGraph parentNode : sub.getParentNodes()) {
      NodeGraph childNode = sub.getChildNode(parentNode);
      assertNotNull(childNode);
      assertSame(parentNode, sub.getParentNode(childNode));
      assertEquals(parentNode.getID(), childNode.getID());
    }
  }

  @Test
  void listOverloadsMapWholeCollections() {
    SubGraph sub = new SubGraph(block);

    List<EdgeGraph> childEdges = sub.getChildEdges(block);
    assertEquals(block.size(), childEdges.size());
    assertEquals(block.size(), sub.getParentEdges(childEdges).size());
    assertTrue(sub.getParentEdges(childEdges).containsAll(block));

    List<NodeGraph> parentNodes = sub.getParentNodes();
    List<NodeGraph> childNodes = sub.getChildNodes(parentNodes);
    assertEquals(parentNodes.size(), childNodes.size());
    assertTrue(sub.getParentNodes(childNodes).containsAll(parentNodes));
  }

  @Test
  @DisplayName("a subgraph is a graph: it answers adjacency questions on its own edges")
  void subGraphAnswersAdjacencyOnItsOwnEdges() {
    SubGraph sub = new SubGraph(block);

    NodeGraph childSouthWest = sub.getChildNode(Fixtures.nodeAt(parent, 0.0, 0.0));
    NodeGraph childSouthEast = sub.getChildNode(Fixtures.nodeAt(parent, 100.0, 0.0));

    assertNotNull(sub.getEdgeBetween(childSouthWest, childSouthEast));
    assertSame(sub.getEdgeBetween(childSouthWest, childSouthEast),
        sub.getEdgeBetween(childSouthEast, childSouthWest));
    // Inside the block every junction has exactly two of the four segments.
    assertEquals(2, childSouthWest.getAdjacentNodes().size());
  }

  @Test
  @DisplayName("the subgraph centrality map is read off the parent nodes")
  void centralityMapComesFromTheParentNodes() {
    List<NodeGraph> parentNodes = parent.getNodes();
    for (int i = 0; i < parentNodes.size(); i++) {
      parentNodes.get(i).setCentrality(i);
    }

    SubGraph sub = new SubGraph(block);
    sub.generateSubGraphCentralityMap();

    Map<NodeGraph, Double> salient = sub.getSubGraphSalientNodes(0.5);
    assertNotNull(salient);
    // The map is keyed on parent nodes, so callers can go back to the parent graph.
    assertTrue(parentNodes.containsAll(salient.keySet()));
  }

  @Test
  @DisplayName("the subgraph salient lookup survives a percentile of 1 and an empty map")
  void subGraphSalientNodesHandleTheEdges() {
    List<NodeGraph> parentNodes = parent.getNodes();
    for (int i = 0; i < parentNodes.size(); i++) {
      parentNodes.get(i).setCentrality(i);
    }

    SubGraph sub = new SubGraph(block);
    sub.generateSubGraphCentralityMap();
    assertEquals(1, sub.getSubGraphSalientNodes(1.0).size());

    assertNull(new SubGraph().getSubGraphSalientNodes(0.5));
  }

  @Test
  void anEmptySubGraphHasNoComponents() {
    SubGraph sub = new SubGraph();
    assertTrue(sub.getEdges().isEmpty());
    assertTrue(sub.getNodes().isEmpty());
  }
}
