
[![Maven Central](https://maven-badges.herokuapp.com/maven-central/uk.ac.liv.gdsl/GeoMason-light/badge.svg?style=plastic)](https://maven-badges.herokuapp.com/maven-central/uk.ac.liv.gdsl/GeoMason-light)
[![GitHub CI](https://github.com/g-filomena/GeoMason-light/actions/workflows/build.yaml/badge.svg)](https://github.com/g-filomena/GeoMason-light/actions/workflows/build.yaml)
[![License](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0.en.html)

## About
GeoMason-light is an optional extension to MASON that adds support for vector and raster geospatial data.

You can find the original manual of Geomason on the [Geomason Site](https://cs.gmu.edu/~eclab/projects/mason/extensions/geomason/).

### Why GeoMason-light?
#### Content:
* It includes a `graph` package that supports a way more advanced manipulation of graph components. This allows the users to develop more complex street network-based urban mobility/traffic agent-based models.
* It includes a `VectorLayer` class that allows for more geometric and spatial operations, usually common for GIS VectorLayers.

#### Support and open software:
* The structure of the repository allows for smoother devleopments and contributions from users, compared both to the MASON and GeoMason source repository.
* It is fully mavenised and available on the [Maven Repository](https://mvnrepository.com).
* It comes with javadoc and sources.

#### Light:
* It comes without all the heavy examples included in the original geomason.
* All the code now refers to the last version of Java Topology Suite (JTS).

### Dependencies:
* mason.
* javatuples.
* jts.

GeoMason-light requires that the jar files for the above dependencies are in the user's java class path or amongst the external libraries of the user's Project in Eclipse.
