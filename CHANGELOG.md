# Changelog

## Version 2.1.0

* [Enhancement] Added `GeoPackageExporter`: writes a `VectorLayer` to a single-file OGC GeoPackage (`.gpkg`), the write counterpart of `GeoPackageImporter`. Uses typed columns and imposes no field-name or field-length limits.
* [Enhancement] Added `GeoJSONExporter`: writes a `VectorLayer` to a GeoJSON `FeatureCollection`, with no external JSON dependency. `toFeatureCollection(layer[, includeProperties])` returns the same document as a string.
* [Enhancement] Added `VectorLayer.writeGPKG(...)` and `VectorLayer.writeGeoJSON(...)` convenience methods, mirroring the existing `readGPKG(...)`.
* [Breaking] Removed `ShapeFileExporter`. Vector features are now written in the open GeoPackage and GeoJSON formats; `ShapeFileImporter` is retained, so existing shapefiles can still be read.
* [Enhancement] Added non-copying `VectorLayer` accessors — `isEmpty()`, `isPopulated()`, `size()` and `geometriesView()` (an unmodifiable view) — so callers that only read the geometries avoid the defensive copy made by `getGeometries()`.
* [Performance] `VectorLayer.getGeometriesFromIDs(Set)` now resolves through a lazily-built id index, so its cost scales with the size of the requested set rather than with the number of geometries in the layer.

## Version 2.0.0

* [Breaking] Started a new versioning line after the legacy `1.x` numbering scheme.
* [Breaking] Updated the project Java baseline from Java 8 to Java 11 to align with current GeoPackage dependencies.
* [Breaking] Declared the `sim:mason:21` dependency as `provided`, since MASON is not available from Maven Central and must be supplied by downstream projects.
* [Fix] Updated Maven metadata, including SCM links and project encoding configuration.
* [Fix] Normalised README dependency examples to use the current release version.
* [Fix] Clarified installation instructions for Maven, Eclipse, local builds, and the external MASON dependency.
* [Enhancement] Added package-level Javadocs for the main API packages:

  * `sim`
  * `sim.field.geo`
  * `sim.graph`
  * `sim.io.geo`
  * `sim.portrayal.geo`
  * `sim.routing`
  * `sim.util.geo`
* [Enhancement] Added a Javadoc style guide under `docs/`.
* [Enhancement] Added GitHub Actions validation for Javadoc generation.
* [Maintenance] Updated GitHub Actions configuration to use Java 11.
* [Maintenance] Improved repository hygiene by excluding local Maven settings, signing material, deployment helpers, generated build output, and local IDE files from version control.

## Version 1.1.7, 1.1.8, 1.1.9
Small fixes and new versioning.

## Version 1.16

* Last release under the legacy `1.x` numbering scheme.
* No detailed changelog entry was recorded for this release.

## Version 1.15

* [Enhancement] Allows loading layers from GeoPackage files.
* [Enhancement] Added A* class.
* [Enhancement] Added further utilities in `GraphUtils.java`.

## Version 1.14

* [Enhancement] Cleaned code; removed or simplified redundant functions.
* [Enhancement] Simplified the `SubGraph` class and added utility methods.
* [Enhancement] Replaced concrete `ArrayList` and `HashMap` declarations with generic `List` and `Map` interfaces where appropriate.

## Version 1.13

* No changelog entry was recorded for this release.

## Version 1.12

* [Fix] Replaced `Bag` objects with `ArrayList` objects, usually for `MasonGeometry` objects.
* [Fix] Renamed and consolidated core classes:

  * `GeomVectorField` merged into `VectorLayer`, previously a subclass of `GeomVectorField`.
  * `GeomGridField` renamed to `GridLayer`.
  * `Field` renamed to `Layer`.
* [Enhancement] Cleaned code inherited from GeoMason.
* [Enhancement] Replaced inefficient loops inherited from GeoMason.
* [Enhancement] Simplified partly redundant functions.
* [Enhancement] Removed specific attributes in the `sim.graph` package for better generalisation.

## Version 1.11

* First release.
