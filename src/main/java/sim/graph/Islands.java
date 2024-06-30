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

import pedSim.routeChoice.Route;

public class Islands {

	Set<NodeGraph> visitedNodes = new HashSet<>();
	private Graph graph;
	List<Set<NodeGraph>> islands = new ArrayList<>();

	public Islands(Graph graph) {
		this.graph = graph;
	}

	public List<Set<NodeGraph>> findDisconnectedIslands(Set<EdgeGraph> edges) {

		Set<NodeGraph> nodes = new HashSet<>(GraphUtils.nodesFromEdges(edges));
		this.visitedNodes = new HashSet<>();
		islands.clear();
		for (NodeGraph currentNode : nodes) {
			if (!visitedNodes.contains(currentNode)) {
				Set<NodeGraph> currentIsland = new HashSet<>();
				dfs(currentNode, currentIsland, nodes, edges);
				islands.add(currentIsland);
			}
		}
		return islands;
	}

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
				Route route = Astar.astarRoute(closestNodes.getValue0(), closestNodes.getValue1(), null);
				bridges.addAll(route.edgesSequence);
			}
			edges.addAll(bridges);
			islands = findDisconnectedIslands(edges);
		}
		return edges;
	}

	private Map.Entry<Pair<NodeGraph, NodeGraph>, Set<NodeGraph>> findClosestPairAcrossAllIslands() {

		return islands.stream()
				.flatMap(island -> islands.stream().filter(otherIsland -> island != otherIsland)
						.flatMap(otherIsland -> island.stream().flatMap(node -> otherIsland.stream().map(
								otherNode -> new AbstractMap.SimpleEntry<>(new Pair<>(node, otherNode), otherIsland)))))
				.min(Comparator.comparingDouble(
						entry -> GraphUtils.nodesDistance(entry.getKey().getValue0(), entry.getKey().getValue1())))
				.orElse(null);
	}

	private static EdgeGraph findConnectingBridge(List<Set<NodeGraph>> islands, Graph graph) {
		return islands.stream()
				.flatMap(island -> islands.stream().filter(otherIsland -> otherIsland != island)
						.flatMap(otherIsland -> island.stream().flatMap(
								node -> otherIsland.stream().map(otherNode -> graph.getEdgeBetween(node, otherNode)))))
				.filter(Objects::nonNull).findAny().orElse(null);
	}

//	public static Set<EdgeGraph> findNecessaryEdges(Set<EdgeGraph> edges, Set<NodeGraph> nodesToPreserve) {
//
//		Set<NodeGraph> allNodes = new HashSet<>(
//				edges.stream().flatMap(edge -> edge.getNodes().stream()).collect(Collectors.toSet()));
//		Set<EdgeGraph> newEdges = new HashSet<>(
//				allNodes.stream().flatMap(node -> node.getEdges().stream()).collect(Collectors.toSet()));
//
//		Map<NodeGraph, Integer> nodeDegrees = new HashMap<>();
//		newEdges.stream().flatMap(edge -> edge.getNodes().stream())
//				.forEach(node -> nodeDegrees.put(node, nodeDegrees.getOrDefault(node, 0) + 1));
//
//		newEdges.stream().filter(edge -> {
//			Set<NodeGraph> nodes = new HashSet<>(edge.getNodes());
//			return nodes.stream()
//					.noneMatch(node -> nodeDegrees.getOrDefault(node, 0) == 1 && !nodesToPreserve.contains(node));
//		}).collect(Collectors.toList());
//
//		return newEdges;
//	}

}
