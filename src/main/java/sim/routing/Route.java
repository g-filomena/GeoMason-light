package sim.routing;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

import org.locationtech.jts.geom.Coordinate;
import org.locationtech.jts.geom.GeometryFactory;
import org.locationtech.jts.geom.LineString;
import org.locationtech.jts.planargraph.DirectedEdge;

import sim.graph.EdgeGraph;
import sim.graph.Graph;
import sim.graph.NodeGraph;

/**
 * A class for storing the sequence of GeomPlanarGraphDirectedEdge in a path and
 * the sequence of NodeWrappers. It supports shortest-path algorithms and
 * provides some utilities.
 */
public class Route {

	// always primal
	public NodeGraph originNode;
	public NodeGraph destinationNode;
	public List<DirectedEdge> directedEdgesSequence = new ArrayList<>();
	public List<NodeGraph> nodesSequence = new ArrayList<>();
	public List<EdgeGraph> edgesSequence = new ArrayList<>();
	public List<NodeGraph> dualNodesSequence = new ArrayList<>();
	public LineString lineString;

	public Map<String, Object> attributes = new HashMap<>();

	private final static GeometryFactory GEOMETRY_FACTORY = new GeometryFactory();

	public Set<NodeGraph> visitedLocations = new HashSet<>();
	public boolean social = false;

	public Route() {

	}

	/**
	 * Constructs a Route object with a given sequence of directed edges.
	 *
	 * @param directedEdgeSequence A list of DirectedEdge objects representing the
	 *                             path.
	 */
	public Route(List<DirectedEdge> directedEdgeSequence) {
		this.directedEdgesSequence = directedEdgeSequence;
	}

	/**
	 * Computes the sequences of nodes and edges for the route, setting the origin
	 * and destination nodes, and obtaining the route geometry.
	 */
	public void computeRouteSequences() {
		nodesSequence();
		edgeSequence();
		originNode = nodesSequence.get(0);
		destinationNode = nodesSequence.get(nodesSequence.size() - 1);
		obtainRouteLineGeometry();

	}

	/**
	 * Returns all the primal nodes traversed in a path.
	 *
	 * @param directedEdgesSequence A sequence of GeomPlanarGraphDirectedEdge
	 *                              representing the path.
	 * @return A list of primal nodes.
	 */
	private void nodesSequence() {
		if (directedEdgesSequence.isEmpty())
			return;

		nodesSequence = directedEdgesSequence.stream().map(directedEdge -> (NodeGraph) directedEdge.getFromNode())
				.collect(Collectors.toList());

		DirectedEdge lastEdge = directedEdgesSequence.get(directedEdgesSequence.size() - 1);
		nodesSequence.add((NodeGraph) lastEdge.getToNode());
	}

	/**
	 * Populates the edgesSequence list with edges from the directed edges sequence.
	 */
	private void edgeSequence() {

		for (DirectedEdge directedEdge : directedEdgesSequence)
			edgesSequence.add((EdgeGraph) directedEdge.getEdge());
	}

	/**
	 * Computes the sequence of dual nodes for the route using the provided dual
	 * network graph.
	 *
	 * @param dualNetwork The graph representing the dual network.
	 */
	public void dualNodesSequence(Graph dualNetwork) {

		directedEdgesSequence.forEach(directedEdge -> {
			EdgeGraph edge = (EdgeGraph) directedEdge.getEdge();
			edgesSequence.add(edge);
			NodeGraph dualNode = dualNetwork.findNode(edge.getDualNode());
			dualNodesSequence.add(dualNode);
		});
	}

	/**
	 * Obtains the geometry of the route line by extracting coordinates from the
	 * edge sequences.
	 */
	private void obtainRouteLineGeometry() {
		List<Coordinate> allCoordinates = new ArrayList<>();

		// Extract coordinates from each LineString segment and add them to the list
		edgesSequence.forEach(edge -> {
			Coordinate fromCoords = edge.getFromNode().getCoordinate();
			List<Coordinate> coordinates = new ArrayList<>(Arrays.asList(edge.getLine().getCoordinates()));

			if (!coordinates.get(0).equals(fromCoords))
				Collections.reverse(coordinates);

			coordinates.stream().filter(coords -> !allCoordinates.contains(coords)).forEach(allCoordinates::add);
		});

		Coordinate[] coordinateArray = allCoordinates.toArray(new Coordinate[0]);
		lineString = GEOMETRY_FACTORY.createLineString(coordinateArray);
	}
}
