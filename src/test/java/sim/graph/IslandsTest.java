package sim.graph;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTimeoutPreemptively;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.time.Duration;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import org.junit.jupiter.api.DisplayName;
import org.junit.jupiter.api.Test;
import sim.testing.Fixtures;

class IslandsTest {

  /** The chain 0-100-200-300-400, with the two middle segments left out of the edge set. */
  private static Set<EdgeGraph> firstAndLastSegment(Graph graph) {
    List<EdgeGraph> edges = new ArrayList<>(graph.getEdges());
    Set<EdgeGraph> subset = new HashSet<>();
    subset.add(edges.get(0));
    subset.add(edges.get(edges.size() - 1));
    return subset;
  }

  @Test
  @DisplayName("a fully connected edge set is a single island")
  void connectedEdgeSetGivesOneIsland() {
    Graph graph = Fixtures.path(5, 100.0);
    List<Set<NodeGraph>> islands =
        new Islands(graph).findDisconnectedIslands(new HashSet<>(graph.getEdges()));

    assertEquals(1, islands.size());
    assertEquals(5, islands.get(0).size());
  }

  @Test
  @DisplayName("omitting the middle segments splits the chain in two")
  void gapsInTheEdgeSetProduceSeveralIslands() {
    Graph graph = Fixtures.path(5, 100.0);
    List<Set<NodeGraph>> islands =
        new Islands(graph).findDisconnectedIslands(firstAndLastSegment(graph));

    assertEquals(2, islands.size());
    for (Set<NodeGraph> island : islands) {
      assertEquals(2, island.size());
    }
  }

  @Test
  @DisplayName("mergeConnectedIslands() bridges islands that the parent graph can join")
  void mergeBridgesIslandsThroughTheParentGraph() {
    Graph graph = Fixtures.path(5, 100.0);
    Set<EdgeGraph> edges = firstAndLastSegment(graph);

    Set<EdgeGraph> merged = new Islands(graph).mergeConnectedIslands(edges);

    assertEquals(graph.getEdges().size(), merged.size());
    assertEquals(1, new Islands(graph).findDisconnectedIslands(merged).size());
  }

  @Test
  @DisplayName("a merge that cannot succeed stops instead of spinning")
  void mergeGivesUpOnUnreachableIslands() {
    // Two chains with no edge and no path between them. Every pass finds no bridging edge and no
    // route across, so nothing is added and the island count never falls: the loop has to give up
    // on its own rather than keep asking the same question.
    Graph graph = Fixtures.twoDisjointPaths();
    Set<EdgeGraph> edges = new HashSet<>(graph.getEdges());
    long before = Islands.incompleteMerges();

    Set<EdgeGraph> merged = assertTimeoutPreemptively(Duration.ofSeconds(10),
        () -> new Islands(graph).mergeConnectedIslands(edges));

    assertEquals(graph.getEdges().size(), merged.size());
    assertEquals(2, new Islands(graph).findDisconnectedIslands(merged).size());
    assertTrue(Islands.incompleteMerges() > before,
        "giving up has to be counted, not silent");
  }

  @Test
  @DisplayName("an already connected edge set is returned untouched")
  void mergeIsANoOpOnAConnectedEdgeSet() {
    Graph graph = Fixtures.grid(3, 3, 100.0);
    Set<EdgeGraph> edges = new HashSet<>(graph.getEdges());
    long before = Islands.incompleteMerges();

    assertEquals(edges.size(), new Islands(graph).mergeConnectedIslands(edges).size());
    assertEquals(before, Islands.incompleteMerges());
  }

  @Test
  @DisplayName("findDisconnectedIslands() can be reused on the same instance")
  void findDisconnectedIslandsResetsBetweenCalls() {
    Graph graph = Fixtures.path(5, 100.0);
    Islands islands = new Islands(graph);

    assertEquals(2, islands.findDisconnectedIslands(firstAndLastSegment(graph)).size());
    assertEquals(1, islands.findDisconnectedIslands(new HashSet<>(graph.getEdges())).size());
  }
}
