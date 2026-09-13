package sim.routing;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertNull;
import static org.junit.jupiter.api.Assertions.assertSame;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.DisplayName;
import org.junit.jupiter.api.Test;
import org.locationtech.jts.planargraph.DirectedEdge;
import sim.graph.EdgeGraph;
import sim.graph.Graph;
import sim.graph.NodeGraph;
import sim.testing.Fixtures;
import sim.util.geo.MasonGeometry;

class RoutingUtilsTest {

  private Graph graph;
  private List<DirectedEdge> chain;

  /** Gives every segment of the graph the dual node that stands for it. */
  private void attachDualNodes() {
    for (EdgeGraph edge : graph.getEdges()) {
      NodeGraph centroid = new NodeGraph(edge.getCoordsCentroid());
      centroid.setMasonGeometry(
          new MasonGeometry(Fixtures.FACTORY.createPoint(edge.getCoordsCentroid())));
      centroid.setPrimalEdge(edge);
      edge.setDualNode(centroid);
    }
  }

  @BeforeEach
  void buildChain() {
    graph = Fixtures.path(4, 100.0);
    attachDualNodes();

    chain = new ArrayList<>();
    NodeGraph previous = Fixtures.nodeAt(graph, 0.0, 0.0);
    for (double x = 100.0; x <= 300.0; x += 100.0) {
      NodeGraph next = Fixtures.nodeAt(graph, x, 0.0);
      chain.add(graph.getDirectedEdgeBetween(previous, next));
      previous = next;
    }
  }

  @Test
  void nodesFromDirectedEdgesSequenceWalksTheWholeChain() {
    List<NodeGraph> nodes = RoutingUtils.getNodesFromDirectedEdgesSequence(chain);

    assertEquals(4, nodes.size());
    assertSame(Fixtures.nodeAt(graph, 0.0, 0.0), nodes.get(0));
    assertSame(Fixtures.nodeAt(graph, 300.0, 0.0), nodes.get(3));
  }

  @Test
  void nodesFromAnEmptySequenceIsEmpty() {
    assertEquals(0,
        RoutingUtils.getNodesFromDirectedEdgesSequence(new ArrayList<DirectedEdge>()).size());
  }

  @Test
  void centroidsFromEdgesSequenceReturnsOneDualNodePerSegment() {
    List<NodeGraph> centroids = RoutingUtils.getCentroidsFromEdgesSequence(chain);

    assertEquals(chain.size(), centroids.size());
    for (int i = 0; i < chain.size(); i++) {
      assertSame(((EdgeGraph) chain.get(i).getEdge()).getDualNode(), centroids.get(i));
    }
  }

  @Test
  @DisplayName("getPrimalJunction() finds the junction two segments share")
  void primalJunctionIsTheSharedNode() {
    NodeGraph shared = Fixtures.nodeAt(graph, 100.0, 0.0);
    EdgeGraph west = graph.getEdgeBetween(Fixtures.nodeAt(graph, 0.0, 0.0), shared);
    EdgeGraph east = graph.getEdgeBetween(shared, Fixtures.nodeAt(graph, 200.0, 0.0));

    assertSame(shared, RoutingUtils.getPrimalJunction(west.getDualNode(), east.getDualNode()));
  }

  @Test
  @DisplayName("getPrimalJunction() is null when the segments do not touch")
  void primalJunctionIsNullForDisjointSegments() {
    EdgeGraph first =
        graph.getEdgeBetween(Fixtures.nodeAt(graph, 0.0, 0.0), Fixtures.nodeAt(graph, 100.0, 0.0));
    EdgeGraph last = graph.getEdgeBetween(Fixtures.nodeAt(graph, 200.0, 0.0),
        Fixtures.nodeAt(graph, 300.0, 0.0));

    assertNull(RoutingUtils.getPrimalJunction(first.getDualNode(), last.getDualNode()));
  }

  @Test
  @DisplayName("a one-segment sequence reports the junction it departs from")
  void previousJunctionOfASingleSegmentIsItsOrigin() {
    assertSame(Fixtures.nodeAt(graph, 0.0, 0.0),
        RoutingUtils.getPreviousJunction(Collections.singletonList(chain.get(0))));
  }

  @Test
  @DisplayName("a longer sequence reports the junction between its last two segments")
  void previousJunctionOfALongerSequenceIsTheLastSharedNode() {
    assertSame(Fixtures.nodeAt(graph, 200.0, 0.0), RoutingUtils.getPreviousJunction(chain));
  }
}
