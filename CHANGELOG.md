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