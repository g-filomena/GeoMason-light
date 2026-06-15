/**
 * Routing support for graph-backed geospatial movement.
 *
 * <p>This package contains shortest-path and route representation utilities built around
 * {@link sim.graph.Graph}, {@link sim.graph.NodeGraph}, and {@link sim.graph.EdgeGraph}. The primary
 * implementation is {@link sim.routing.Astar}, supported by route containers and node wrappers.</p>
 *
 * <h2>Routing model</h2>
 * <ul>
 *   <li>{@link sim.routing.Astar}: A* path search over graph nodes and edges.</li>
 *   <li>{@link sim.routing.Route}: ordered route result object for simulation agents.</li>
 *   <li>{@link sim.routing.NodeWrapper}: search-state wrapper around graph nodes.</li>
 *   <li>{@link sim.routing.RoutingUtils}: shared helper methods for route construction.</li>
 * </ul>
 *
 * <p>Routing costs should be documented by callers when they differ from geometric edge length,
 * because route choice and behavioural models depend on the meaning of those weights.</p>
 *
 * @see sim.graph.Graph
 * @see sim.graph.EdgeGraph
 * @see sim.graph.NodeGraph
 */
package sim.routing;
