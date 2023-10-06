## About
Geomason-light is an optional extension to MASON that adds support for vector and raster geospatial data.

You can find the original Manual on the [Geomason Site](https://cs.gmu.edu/~eclab/projects/mason/extensions/geomason/)

### Why geomason-light?
#### Content:
* It includes a graph package that supports a way more advanced manipulation of graph components. This allows the users to developed more complex street network-based urban mobility/traffic simulation agent-based models.

#### Support and open software:
* The structure of the repository allow for smoother devleopments and contributions from users (the mason and geomason repository appear rather close at the moment).
* It is fully mavenised and available on the [Maven Repository](https://mvnrepository.com).
* It comes with javadoc and sources.

#### Light:
* It comes without all the heavy examples included in the original geomason.
* All the code now refers to the last version of Java Topology Suite (JTS).

### Dependencies:
* mason;
* javatuples;
* jts.

geomason-light requires that the jar files for the above dependencies are in the user's java class path or amongst the external libraries of the user's Project in Eclipse.
