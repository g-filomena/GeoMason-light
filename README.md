

[![Maven Central](https://maven-badges.herokuapp.com/maven-central/uk.ac.liv.gdsl/GeoMason-light/badge.svg?style=plastic)](https://maven-badges.herokuapp.com/maven-central/uk.ac.liv.gdsl/GeoMason-light)
[![GitHub CI](https://github.com/g-filomena/GeoMason-light/actions/workflows/build.yaml/badge.svg)](https://github.com/g-filomena/GeoMason-light/actions/workflows/build.yaml)
[![License](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0.en.html)

## About
`GeoMason-light` is an optional extension to MASON that adds support for vector and raster geospatial data.

You can find the original manual of `Geomason` on the [Geomason Site](https://cs.gmu.edu/~eclab/projects/mason/extensions/geomason/).

### Why GeoMason-light?
#### Content:
* It includes a `graph` package that supports a way more advanced manipulation of graph components. This allows the users to develop flexible street network-based urban mobility/traffic agent-based models.
* It explicitly allows the creation of `SubGraphs`.
* It includes a `VectorLayer` class that allows for more geometric and spatial operations, usually common for GIS VectorLayers.
* It supports loading GeoPackages through the library [Geopackage](https://github.com/ngageoint/geopackage-java).

#### Open software:
* The structure of the repository allows for smoother devleopments and contributions from users, compared both to the MASON and GeoMason source repositories.
* It is fully mavenised and available on the [Maven Repository](https://mvnrepository.com).
* It comes with javadoc and sources.

#### Light:
* It comes without all the heavy examples included in the original geomason.
* All the code now refers to the last version of Java Topology Suite (JTS).

### Main Dependencies:
* `mason`
* `javatuples`
* `jts`

# Installation Instructions

Follow the steps below to install the necessary tools, and run the GeoMason-light project locally. Normally, packages relying on GeoMason-light are already configured to install GeoMason and its dependencies automatically.

## Option 1: Installing GeoMason-light in Eclipse

1.  Ensure that **Maven** is installed and configured. If not, you can install it from the **Eclipse Marketplace**. 
2. Create a New Maven Project (If you don't already have one).  In Eclipse, go to **File** → **New** → **Maven Project**.  Follow the prompts to create a new Maven project.
3. Add GeoMason-light Dependency to the `pom.xml` file of your project:
```
<dependency> <groupId>uk.ac.liv.gdsl</groupId> 
<artifactId>GeoMason-light</artifactId> 
<version>1.17</version> </dependency>
```
Alternatively, if you are not using Maven, add `GeoMason-light-1.1.7` to your Java Project's **Class Path** (right-click on project, **Build Path...**, **Configure Path**) through the **Add External JARs...**

You can now start using GeoMason-light in your project.     
   
### Option 2: Installing GeoMason-light into Maven (Central Repository)
GeoMason-light is available in **Maven Central**, so you can directly install it into your local Maven repository. This allows you to use it in other Maven projects, including those not in Eclipse.

1. Download Maven from the [official Maven website](https://maven.apache.org/download.cgi) and extract Maven to a directory on your machine.
3.  Set up Environment Variables: Add Maven's `bin` directory to your system's `PATH` environment variable
4. Add GeoMason-light Dependency to Your `pom.xml`, inside the `<dependencies>` section:

```
<dependency> <groupId>uk.ac.liv.gdsl</groupId> 
<artifactId>GeoMason-light</artifactId> 
<version>1.17</version> </dependency>
```

Since `GeoMason-light` is available in Maven Central, Maven will automatically resolve and download the dependency for you when you run a build. To install the dependency:

1.  Open a terminal or command prompt.
2.  Navigate to the root directory of your Maven project.
3.  Run the following Maven command to download the dependencies and build your project:
   
`mvn clean install` 

4. Verify the Dependency. After building your project, check that **GeoMason-light** has been added to the local Maven repository (typically located at `C:\Users\<user>\.m2\repository` on Windows). Once the dependency is added to your `pom.xml`, you can start using GeoMason-light in your Java classes:

