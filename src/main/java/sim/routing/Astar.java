package sim.routing;

import java.util.ArrayList;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.PriorityQueue;
import java.util.Set;

import org.locationtech.jts.geom.Coordinate;
import org.locationtech.jts.planargraph.DirectedEdge;

import sim.graph.EdgeGraph;
import sim.graph.Graph;
import sim.graph.NodeGraph;
import sim.util.geo.GeometryUtilities;

public class Astar {

	/**
	 * Finds the shortest route between two nodes using the A* algorithm, avoiding
	 * specified edges.
	 *
	 * @param originNode      the starting node for the route
	 * @param destinationNode the destination node for the route
	 * @param graph           the graph on which to calculate the route
	 * @param edgesToAvoid    a set of edges to avoid during route calculation
	 * @return a Route object representing the shortest path from originNode to
	 *         destinationNode, or null if no path is found
	 */
	public static Route astarRoute(NodeGraph originNode, NodeGraph destinationNode, Graph graph,
			Set<EdgeGraph> edgesToAvoid) {

		if (edgesToAvoid == null)
			edgesToAvoid = new HashSet<>();

		// Data structures for A* algorithm
		Map<NodeGraph, NodeWrapper> nodeWrappersMap = new HashMap<>();
		PriorityQueue<NodeWrapper> openSet = new PriorityQueue<>(Comparator.comparingDouble(n -> n.fx));
		Set<NodeGraph> closedSet = new HashSet<>();

		// Initialize the origin node wrapper
		NodeWrapper originNodeWrapper = getNodeWrapper(originNode, nodeWrappersMap);
		originNodeWrapper.gx = 0.0;
		originNodeWrapper.hx = heuristic(originNode, destinationNode);
		originNodeWrapper.fx = originNodeWrapper.hx;
		openSet.add(originNodeWrapper);

		while (!openSet.isEmpty()) {

			NodeWrapper currentNodeWrapper = openSet.poll();
			NodeGraph currentNode = currentNodeWrapper.node;

			if (currentNode.equals(destinationNode))
				return reconstructPath(currentNodeWrapper, originNode, destinationNode);
			closedSet.add(currentNode);

			for (EdgeGraph edge : currentNode.getEdges()) {

				if (edgesToAvoid.contains(edge))
					continue;
				NodeGraph otherNode = edge.getOtherNode(currentNode);
				if (closedSet.contains(otherNode))
					continue;

				NodeWrapper nextNodeWrapper = getNodeWrapper(otherNode, nodeWrappersMap);
				double tentativeGx = currentNodeWrapper.gx + edge.getLength();

				if (tentativeGx < nextNodeWrapper.gx) {
					nextNodeWrapper.previousWrapper = currentNodeWrapper;
					nextNodeWrapper.directedEdgeFrom = graph.getDirectedEdgeBetween(currentNode, otherNode);
					nextNodeWrapper.gx = tentativeGx;
					nextNodeWrapper.hx = heuristic(otherNode, destinationNode); // Update hx value
					nextNodeWrapper.fx = tentativeGx + nextNodeWrapper.hx;

					// Only add to open set if it's not already there
					if (!openSet.contains(nextNodeWrapper))
						openSet.add(nextNodeWrapper);
					else {
						// Update the position in the priority queue
						openSet.remove(nextNodeWrapper);
						openSet.add(nextNodeWrapper);
					}
				}
			}
		}
		return null;
	}

	/**
	 * Retrieves the NodeWrapper associated with the given node, creating it if it
	 * doesn't exist.
	 *
	 * @param node            the node to get the wrapper for
	 * @param nodeWrappersMap a map storing NodeWrappers for each node
	 * @return the NodeWrapper for the given node
	 */
	private static NodeWrapper getNodeWrapper(NodeGraph node, Map<NodeGraph, NodeWrapper> nodeWrappersMap) {
		return nodeWrappersMap.computeIfAbsent(node, NodeWrapper::new);
	}

	/**
	 * Reconstructs the path from the origin node to the destination node based on
	 * the node wrappers.
	 *
	 * @param nodeWrapper     the NodeWrapper at the destination node
	 * @param originNode      the starting node for the route
	 * @param destinationNode the destination node for the route
	 * @return a Route object representing the path from originNode to
	 *         destinationNode
	 */
	private static Route reconstructPath(NodeWrapper nodeWrapper, NodeGraph originNode, NodeGraph destinationNode) {

		Route route = new Route();
		List<DirectedEdge> directedEdgesSequence = new ArrayList<>();

		NodeWrapper currentNodeWrapper = nodeWrapper;
		while (currentNodeWrapper.previousWrapper != null) {
			directedEdgesSequence.add(0, currentNodeWrapper.directedEdgeFrom);
			currentNodeWrapper = currentNodeWrapper.previousWrapper;
		}
		route.directedEdgesSequence = directedEdgesSequence;
		if (!route.directedEdgesSequence.isEmpty())
			route.computeRouteSequences();
		return route;
	}

	/**
	 * Calculates the heuristic distance between two nodes using Euclidean distance.
	 *
	 * @param node      the first node
	 * @param otherNode the second node
	 * @return the heuristic distance between the two nodes
	 */
	private static double heuristic(NodeGraph node, NodeGraph otherNode) {

		Coordinate nodeCoords = node.getCoordinate();
		Coordinate otherNodeCoords = otherNode.getCoordinate();
		return GeometryUtilities.euclideanDistance(nodeCoords, otherNodeCoords);
	}
}
