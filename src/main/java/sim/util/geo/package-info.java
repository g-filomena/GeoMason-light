/**
 * Utility classes for geometry manipulation, attributes, angles, CSV handling, and coordinate moves.
 *
 * <p>The package contains small reusable helpers used throughout GeoMason-light. The central value
 * object is {@link sim.util.geo.MasonGeometry}, which wraps JTS geometries with simulation-facing
 * attributes.</p>
 *
 * <h2>Utility groups</h2>
 * <ul>
 *   <li>{@link sim.util.geo.GeometryUtilities}: geometry construction and transformation helpers.</li>
 *   <li>{@link sim.util.geo.Angles}: angular calculations used by graph and routing workflows.</li>
 *   <li>{@link sim.util.geo.AttributeValue}: typed attribute container for feature metadata.</li>
 *   <li>{@link sim.util.geo.CSVUtils}: CSV parsing and writing helpers.</li>
 *   <li>{@link sim.util.geo.PointMoveTo}: coordinate-sequence transformation support.</li>
 * </ul>
 *
 * <p>Prefer these helpers over duplicating geometry or attribute logic in higher-level packages.</p>
 *
 * @see sim.field.geo.VectorLayer
 * @see sim.graph.EdgeGraph
 */
package sim.util.geo;
