/**
 * Import and export utilities for geospatial file formats used by GeoMason-light.
 *
 * <p>This package contains readers and writers for raster and vector formats, including ArcInfo ASCII
 * grids, GeoPackages, and GeoJSON. Importers normally populate {@link sim.field.geo.GridLayer}
 * or {@link sim.field.geo.VectorLayer} instances with geometry and attribute data that can then be
 * used by simulations, graphs, and portrayals. Vector export uses the open GeoPackage and GeoJSON
 * formats; ESRI Shapefiles can still be read but are no longer written.</p>
 *
 * <h2>Supported workflows</h2>
 * <ul>
 *   <li>Read raster grids with {@link sim.io.geo.ArcInfoASCGridImporter}.</li>
 *   <li>Write raster grids with {@link sim.io.geo.ArcInfoASCGridExporter}.</li>
 *   <li>Read vector features with {@link sim.io.geo.GeoPackageImporter} and
 *       {@link sim.io.geo.ShapeFileImporter}.</li>
 *   <li>Write vector features with {@link sim.io.geo.GeoPackageExporter} and
 *       {@link sim.io.geo.GeoJSONExporter}.</li>
 * </ul>
 *
 * <p>Callers should validate coordinate reference system assumptions outside these importers when
 * mixing layers from different sources.</p>
 *
 * @see sim.field.geo.VectorLayer
 * @see sim.field.geo.GridLayer
 * @see sim.util.geo.MasonGeometry
 */
package sim.io.geo;
