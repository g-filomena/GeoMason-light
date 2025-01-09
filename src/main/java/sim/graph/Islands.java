package sim.graph;

import java.util.AbstractMap;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Objects;
import java.util.Set;
import java.util.Stack;

import org.javatuples.Pair;

import sim.routing.Astar;
import sim.routing.Route;

public class Islands {

	Set<NodeGraph> visitedNodes = new HashSet<>();
	private Graph graph;
	List<Set<NodeGraph>> islands = new ArrayList<>();

	/**
	 * Constructs an Islands object for finding and merging disconnected islands in
	 * the given graph.
	 *
	 * @param graph the graph to analyze for islands
	 */
	public Islands(Graph graph) {
		this.graph = graph;
	}

	/**
	 * Finds and returns a list of sets, where each set contains nodes that form a
	 * disconnected island in the graph.
	 *
	 * @param edges the set of edges to consider when finding disconnected islands
	 * @return a list of sets, each representing a disconnected island of nodes
	 */
	public List<Set<NodeGraph>> findDisconnectedIslands(Set<EdgeGraph> edges) {
		Set<NodeGraph> nodes = new HashSet<>(GraphUtils.nodesFromEdges(edges));
		this.visitedNodes = new HashSet<>();
		islands.clear();

		nodes.parallelStream().forEach(currentNode -> {
			if (!visitedNodes.contains(currentNode)) {
				synchronized (visitedNodes) { // Ensure thread-safe updates
					if (!visitedNodes.contains(currentNode)) {
						Set<NodeGraph> currentIsland = new HashSet<>();
						dfs(currentNode, currentIsland, nodes, edges);
						synchronized (islands) {
							islands.add(currentIsland);
						}
					}
				}
			}
		});

		return islands;
	}

	/**
	 * Performs a depth-first search (DFS) to explore and mark all nodes connected
	 * to the given node.
	 *
	 * @param node          the starting node for the DFS
	 * @param currentIsland the set to store nodes that are part of the current
	 *                      island
	 * @param nodes         the set of all nodes to consider in the DFS
	 * @param edges         the set of edges to consider in the DFS
	 */
	private void dfs(NodeGraph node, Set<NodeGraph> currentIsland, Set<NodeGraph> nodes, Set<EdgeGraph> edges) {

		Stack<NodeGraph> stack = new Stack<>();
		stack.push(node);

		while (!stack.isEmpty()) {
			NodeGraph currentNode = stack.pop();

			if (!visitedNodes.contains(currentNode)) {
				visitedNodes.add(currentNode);
				currentIsland.add(currentNode);
				List<NodeGraph> adjacentNodes = new ArrayList<>(currentNode.getAdjacentNodes());
				adjacentNodes.retainAll(nodes);
				currentNode.getAdjacentNodes().stream()
						.filter(adjacentNode -> edges.contains(graph.getEdgeBetween(currentNode, adjacentNode)))
						.forEach(stack::push);
			}
		}
	}

	/**
	 * Merges disconnected islands in the graph by adding edges to connect them,
	 * ensuring all nodes are connected.
	 *
	 * @param edges the set of edges to consider when merging islands
	 * @return the updated set of edges with added bridges to connect islands
	 */
	public Set<EdgeGraph> mergeConnectedIslands(Set<EdgeGraph> edges) {

		Set<EdgeGraph> bridges = new HashSet<>();
		islands = findDisconnectedIslands(edges);
		if (islands.size() == 1)
			return edges;

		// Connect the remaining islands purely based on distance
		while (islands.size() > 1) {
			EdgeGraph connectingEdge = findConnectingBridge(islands, graph);
			if (connectingEdge != null)
				bridges.add(connectingEdge);
			else {
				Map.Entry<Pair<NodeGraph, NodeGraph>, Set<NodeGraph>> closestPairWithIsland = findClosestPairAcrossAllIslands();
				Pair<NodeGraph, NodeGraph> closestNodes = closestPairWithIsland.getKey();
				Route route = Astar.astarRoute(closestNodes.getValue0(), closestNodes.getValue1(), graph, null);
				bridges.addAll(route.edgesSequence);
			}
			edges.addAll(bridges);
			islands = findDisconnectedIslands(edges);
		}
		return edges;
	}

	/**
	 * Finds the closest pair of nodes across all islands and returns the pair along
	 * with the island of the second node.
	 *
	 * @return a map entry containing the closest pair of nodes and the island of
	 *         the second node
	 */
	private Map.Entry<Pair<NodeGraph, NodeGraph>, Set<NodeGraph>> findClosestPairAcrossAllIslands() {
		return islands.parallelStream()
				.flatMap(island -> islands.parallelStream().filter(otherIsland -> island != otherIsland)
						.flatMap(otherIsland -> island.stream().flatMap(node -> otherIsland.stream().map(
								otherNode -> new AbstractMap.SimpleEntry<>(new Pair<>(node, otherNode), otherIsland)))))
				.min(Comparator.comparingDouble(
						entry -> GraphUtils.nodesDistance(entry.getKey().getValue0(), entry.getKey().getValue1())))
				.orElse(null);
	}

	/**
	 * Finds an edge that can act as a bridge to connect two disconnected islands in
	 * the graph.
	 *
	 * @param islands the list of sets of nodes representing the disconnected
	 *                islands
	 * @param graph   the graph to analyze for potential connecting bridges
	 * @return an edge that connects two islands, or null if no such edge exists
	 */
	private static EdgeGraph findConnectingBridge(List<Set<NodeGraph>> islands, Graph graph) {
		return islands.stream()
				.flatMap(island -> islands.stream().filter(otherIsland -> otherIsland != island)
						.flatMap(otherIsland -> island.stream().flatMap(
								node -> otherIsland.stream().map(otherNode -> graph.getEdgeBetween(node, otherNode)))))
				.filter(Objects::nonNull).findAny().orElse(null);
	}
}
