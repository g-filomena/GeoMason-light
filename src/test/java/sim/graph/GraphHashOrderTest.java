package sim.graph;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertNotEquals;
import static org.junit.jupiter.api.Assertions.assertNotNull;
import static org.junit.jupiter.api.Assertions.assertNotSame;
import static org.junit.jupiter.api.Assertions.assertSame;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import org.junit.jupiter.api.DisplayName;
import org.junit.jupiter.api.Test;
import org.locationtech.jts.geom.Coordinate;
import sim.testing.Fixtures;

/**
 * The hashing contract that makes a seeded run reproduce on a different machine.
 *
 * <p>{@code NodeGraph} and {@code EdgeGraph} hash on their geometry rather than on identity.
 * HotSpot draws an identity hash from a per-JVM generator, so a {@code HashMap} or {@code HashSet}
 * keyed on one would iterate in an order differing between JVM builds, and any decision taken by
 * walking such a collection in order would reach a different answer on a different machine from
 * the same seed.
 */
class GraphHashOrderTest {

  @Test
  @DisplayName("a node's hash comes from its coordinate, not from its identity")
  void nodeHashIsCoordinateDerived() {
    Graph graph = Fixtures.grid(3, 3, 100.0);
    NodeGraph node = Fixtures.nodeAt(graph, 100.0, 100.0);
    assertNotNull(node);
    assertEquals(new Coordinate(100.0, 100.0).hashCode(), node.hashCode());
  }

  @Test
  @DisplayName("a node's hash survives setID, which is called after the graph is built")
  void nodeHashIsStableAcrossSetID() {
    // The trap this guards: Graph.generateAdjacencyMatrix stores Pair<NodeGraph, NodeGraph> keys
    // while the graph is being built, and callers assign IDs afterwards. A hash derived from the ID
    // would move every one of those keys to a bucket it is not stored in, and adjacency lookups
    // would start returning null.
    Graph graph = Fixtures.grid(3, 3, 100.0);
    NodeGraph node = Fixtures.nodeAt(graph, 100.0, 100.0);

    int before = node.hashCode();
    node.setID(4242);
    assertEquals(before, node.hashCode(), "hash must not depend on a field assigned after building");

    node.setRegionID(7);
    assertEquals(before, node.hashCode(), "hash must not depend on the region either");
  }

  @Test
  @DisplayName("adjacency lookups still resolve after IDs are assigned")
  void adjacencySurvivesSetID() {
    Graph graph = Fixtures.grid(3, 3, 100.0);
    NodeGraph from = Fixtures.nodeAt(graph, 0.0, 0.0);
    NodeGraph to = Fixtures.nodeAt(graph, 100.0, 0.0);
    assertNotNull(graph.getEdgeBetween(from, to));

    int id = 0;
    for (NodeGraph node : graph.getNodes()) {
      node.setID(id++);
    }
    assertNotNull(
        graph.getEdgeBetween(from, to), "adjacency must still resolve once IDs are assigned");
  }

  @Test
  @DisplayName("equals stays identity, so distinct nodes at one coordinate remain distinct")
  void equalsRemainsIdentity() {
    NodeGraph one = new NodeGraph(new Coordinate(5.0, 5.0));
    NodeGraph other = new NodeGraph(new Coordinate(5.0, 5.0));

    assertNotSame(one, other);
    assertNotEquals(one, other, "equals must not become coordinate-based");
    assertEquals(one.hashCode(), other.hashCode(), "colliding hashes are allowed; they share a bucket");

    Set<NodeGraph> set = new HashSet<>();
    set.add(one);
    set.add(other);
    assertEquals(2, set.size(), "a hash collision must not merge two distinct nodes");
  }

  @Test
  @DisplayName("a node's hash does not vary between equivalent graphs built separately")
  void hashIsReproducibleAcrossBuilds() {
    // Two graphs built identically give equal hashes for corresponding nodes. Under an identity
    // hash this held only by accident, and never across JVMs.
    Graph first = Fixtures.grid(3, 3, 100.0);
    Graph second = Fixtures.grid(3, 3, 100.0);

    List<Integer> firstHashes = new ArrayList<>();
    List<Integer> secondHashes = new ArrayList<>();
    for (NodeGraph node : first.getNodes()) {
      firstHashes.add(node.hashCode());
    }
    for (NodeGraph node : second.getNodes()) {
      secondHashes.add(node.hashCode());
    }
    assertEquals(firstHashes, secondHashes);
  }

  @Test
  @DisplayName("an edge's hash comes from its centroid and survives setID")
  void edgeHashIsCentroidDerivedAndStable() {
    Graph graph = Fixtures.grid(3, 3, 100.0);
    EdgeGraph edge = graph.getEdges().get(0);

    int before = edge.hashCode();
    edge.setID(99);
    assertEquals(before, edge.hashCode());
    assertEquals(edge.getCoordsCentroid().hashCode(), edge.hashCode());
  }

  @Test
  @DisplayName("a set of graph objects iterates in the same order for identical graphs")
  void iterationOrderMatchesForIdenticalGraphs() {
    Graph first = Fixtures.grid(4, 4, 50.0);
    Graph second = Fixtures.grid(4, 4, 50.0);

    Set<NodeGraph> firstSet = new HashSet<>(first.getNodes());
    Set<NodeGraph> secondSet = new HashSet<>(second.getNodes());

    List<Coordinate> firstOrder = new ArrayList<>();
    List<Coordinate> secondOrder = new ArrayList<>();
    for (NodeGraph node : firstSet) {
      firstOrder.add(node.getCoordinate());
    }
    for (NodeGraph node : secondSet) {
      secondOrder.add(node.getCoordinate());
    }
    // This is the property the whole change exists for: the walk order of a HashSet of graph
    // objects is now a function of the graph, not of the JVM that happens to be running.
    assertEquals(firstOrder, secondOrder);
  }
}
