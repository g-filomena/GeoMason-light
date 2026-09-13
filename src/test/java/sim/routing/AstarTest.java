package sim.routing;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertNotNull;
import static org.junit.jupiter.api.Assertions.assertNull;
import static org.junit.jupiter.api.Assertions.assertSame;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.util.Collections;
import java.util.HashSet;
import java.util.Set;
import org.junit.jupiter.api.DisplayName;
import org.junit.jupiter.api.Test;
import sim.graph.EdgeGraph;
import sim.graph.Graph;
import sim.graph.NodeGraph;
import sim.testing.Fixtures;

class AstarTest {

  @Test
  @DisplayName("on a grid the route is as long as the Manhattan distance")
  void routesAcrossAGridAreShortest() {
    Graph graph = Fixtures.grid(3, 3, 100.0);
    NodeGraph origin = Fixtures.nodeAt(graph, 0.0, 0.0);
    NodeGraph destination = Fixtures.nodeAt(graph, 200.0, 200.0);

    Route route = new Astar().astarRoute(origin, destination, graph, null);

    assertNotNull(route);
    assertEquals(4, route.edgesSequence.size());
    assertEquals(5, route.nodesSequence.size());
    assertEquals(400.0, route.getLength(), 1e-9);
    assertSame(origin, route.originNode);
    assertSame(destination, route.destinationNode);
  }

  @Test
  @DisplayName("along a chain the route is the whole chain")
  void routesAlongAChainFollowEverySegment() {
    Graph graph = Fixtures.path(5, 100.0);
    NodeGraph origin = Fixtures.nodeAt(graph, 0.0, 0.0);
    NodeGraph destination = Fixtures.nodeAt(graph, 400.0, 0.0);

    Route route = new Astar().astarRoute(origin, destination, graph, null);

    assertNotNull(route);
    assertEquals(4, route.edgesSequence.size());
    assertEquals(400.0, route.getLength(), 1e-9);
    assertTrue(route.edgesSequence.containsAll(graph.getEdges()));
  }

  @Test
  @DisplayName("a route to the origin itself is empty rather than null")
  void routeToSelfIsEmpty() {
    Graph graph = Fixtures.grid(3, 3, 100.0);
    NodeGraph origin = Fixtures.nodeAt(graph, 0.0, 0.0);

    Route route = new Astar().astarRoute(origin, origin, graph, null);

    assertNotNull(route);
    assertTrue(route.directedEdgesSequence.isEmpty());
    assertTrue(route.edgesSequence.isEmpty());
    assertEquals(0.0, route.getLength(), 0.0);
  }

  @Test
  @DisplayName("edges to avoid are not traversed")
  void avoidedEdgesAreDetoured() {
    Graph graph = Fixtures.grid(3, 3, 100.0);
    NodeGraph origin = Fixtures.nodeAt(graph, 0.0, 0.0);
    NodeGraph destination = Fixtures.nodeAt(graph, 100.0, 0.0);
    EdgeGraph direct = graph.getEdgeBetween(origin, destination);

    Set<Integer> avoid = new HashSet<>(Collections.singletonList(direct.getID()));
    Route route = new Astar().astarRoute(origin, destination, graph, avoid);

    assertNotNull(route);
    assertTrue(route.edgesSequence.stream().noneMatch(edge -> edge.getID() == direct.getID()));
    // The detour runs round the block: three sides instead of one.
    assertEquals(300.0, route.getLength(), 1e-9);
  }

  @Test
  @DisplayName("a destination cut off by the avoided edges is unreachable")
  void unreachableDestinationReturnsNull() {
    Graph graph = Fixtures.path(3, 100.0);
    NodeGraph origin = Fixtures.nodeAt(graph, 0.0, 0.0);
    NodeGraph middle = Fixtures.nodeAt(graph, 100.0, 0.0);
    NodeGraph destination = Fixtures.nodeAt(graph, 200.0, 0.0);

    Set<Integer> avoid =
        new HashSet<>(Collections.singletonList(graph.getEdgeBetween(middle, destination).getID()));

    assertNull(new Astar().astarRoute(origin, destination, graph, avoid));
  }

  @Test
  @DisplayName("a destination in another component is unreachable")
  void disconnectedDestinationReturnsNull() {
    Graph graph = Fixtures.twoDisjointPaths();
    NodeGraph origin = Fixtures.nodeAt(graph, 0.0, 0.0);
    NodeGraph destination = Fixtures.nodeAt(graph, 10000.0, 0.0);

    assertNull(new Astar().astarRoute(origin, destination, graph, null));
  }

  @Test
  @DisplayName("the same query answers the same way twice")
  void routingIsDeterministic() {
    Graph graph = Fixtures.grid(4, 4, 100.0);
    NodeGraph origin = Fixtures.nodeAt(graph, 0.0, 0.0);
    NodeGraph destination = Fixtures.nodeAt(graph, 300.0, 300.0);

    Route first = new Astar().astarRoute(origin, destination, graph, null);
    Route second = new Astar().astarRoute(origin, destination, graph, null);

    assertEquals(first.edgesSequence, second.edgesSequence);
    assertEquals(first.getLength(), second.getLength(), 0.0);
  }
}
