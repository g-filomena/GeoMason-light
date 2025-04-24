Changelog
=========

Version 1.11
------------
- First release.

Version 1.12
------------
- [Fix] Replaced `Bag` objects with `ArrayList` objects (usually for `MasonGeometry` objects).
- [Fix] Class names: `GeomVectorField` merged into --> `VectorLayer` (previously a subclass of `GeomVectorField`), `GeomGridField` --> `GridLayer`, `Field` --> `Layer`.
- [Enhancement] Cleaned code; replaced inefficient for loops inherited from GeoMason; simplified partly redundant functions.
- [Enhancement] Removed specific attributes in the `sim.graph` package for better generalisation.

Version 1.14
------------
- [Enhancement] Cleaned code; removed/simplified redundant functions. Simplified class SubGraph and added some utilities.
- [Enhancement] Replaced ArrayList/HashMap with generic types List/Map

Version 1.15
------------
- [Enhancement] Allows loading layers from Geopackage files.
- [Enhancement] Added A* Class.
- [Enhancement] Added further utilities in GraphUtils.java