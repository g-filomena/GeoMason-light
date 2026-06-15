/**
 * Geospatial field abstractions for storing raster and vector data in MASON simulations.
 *
 * <p>This package groups the common {@link sim.field.geo.Layer} contract with concrete layer
 * implementations such as {@link sim.field.geo.GridLayer} for raster grids and
 * {@link sim.field.geo.VectorLayer} for JTS-backed feature collections.</p>
 *
 * <h2>Typical responsibilities</h2>
 * <ul>
 *   <li>Store simulation geometries and grid cells in spatially indexed structures.</li>
 *   <li>Keep layer envelopes aligned with the simulation coordinate reference system.</li>
 *   <li>Expose layer data to importers, exporters, graph builders, and portrayals.</li>
 * </ul>
 *
 * <p>Attribute-bearing geometries should usually be represented with
 * {@link sim.util.geo.MasonGeometry} before being added to vector layers.</p>
 *
 * @see sim.util.geo.MasonGeometry
 * @see sim.io.geo.ShapeFileImporter
 * @see sim.io.geo.GeoPackageImporter
 */
package sim.field.geo;
