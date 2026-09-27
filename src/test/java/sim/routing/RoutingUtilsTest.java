package sim.routing;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertNull;
import static org.junit.jupiter.api.Assertions.assertSame;

import java.util.ArrayList;
import java.util.Arrays;
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

  /**
   * A road from (-100,0) to junction A (0,0), then two parallel streets between A and B (100,0): a
   * straight one and a crescent. Parallel streets share both ends, so the junction between them is
   * the one the walk did not arrive by.
   */
  private List<EdgeGraph> approachThenParallelStreets() {
    graph = Fixtures.graphOf(Arrays.asList(Fixtures.segment(-100, 0, 0, 0),
        Fixtures.segment(0, 0, 100, 0), Fixtures.polyline(0, 0, 50, 60, 100, 0)));
    attachDualNodes();
    NodeGraph a = Fixtures.nodeAt(graph, 0.0, 0.0);
    NodeGraph b = Fixtures.nodeAt(graph, 100.0, 0.0);
    List<EdgeGraph> parallel = graph.getEdgesBetween(a, b);
    return Arrays.asList(graph.getEdgeBetween(Fixtures.nodeAt(graph, -100.0, 0.0), a),
        parallel.get(0), parallel.get(1));
  }

  @Test
  @DisplayName("between parallel streets, the junction is the far end from the arrival")
  void primalJunctionOfParallelStreetsIsTheFarEnd() {
    List<EdgeGraph> streets = approachThenParallelStreets();
    NodeGraph a = Fixtures.nodeAt(graph, 0.0, 0.0);
    NodeGraph b = Fixtures.nodeAt(graph, 100.0, 0.0);
    NodeGraph straight = streets.get(1).getDualNode();
    NodeGraph crescent = streets.get(2).getDualNode();

    assertSame(b, RoutingUtils.getPrimalJunction(straight, crescent, a));
    assertSame(a, RoutingUtils.getPrimalJunction(straight, crescent, b));
    // without the arrival, parallel streets fall back to the first segment's from-node
    assertSame(streets.get(1).getFromNode(),
        RoutingUtils.getPrimalJunction(straight, crescent, null));
  }

  @Test
  @DisplayName("the previous junction of a walk onto a parallel street is the far end")
  void previousJunctionWalksPastParallelStreets() {
    List<EdgeGraph> streets = approachThenParallelStreets();
    List<DirectedEdge> walk = new ArrayList<>();
    for (EdgeGraph street : streets) {
      walk.add(street.getDirEdge(0));
    }

    // approach -> A -> straight -> B -> crescent: the last junction crossed is B
    assertSame(Fixtures.nodeAt(graph, 100.0, 0.0), RoutingUtils.getPreviousJunction(walk));
  }
}
