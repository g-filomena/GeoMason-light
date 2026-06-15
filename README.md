# GeoMason-light

[![Maven Central](https://img.shields.io/maven-central/v/uk.ac.liv.gdsl/GeoMason-light.svg)](https://mvnrepository.com/artifact/uk.ac.liv.gdsl/GeoMason-light)
[![GitHub CI](https://github.com/g-filomena/GeoMason-light/actions/workflows/build.yaml/badge.svg)](https://github.com/g-filomena/GeoMason-light/actions/workflows/build.yaml)
[![Javadoc](https://github.com/g-filomena/GeoMason-light/actions/workflows/javadoc.yml/badge.svg)](https://github.com/g-filomena/GeoMason-light/actions/workflows/javadoc.yml)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0.en.html)

## About

`GeoMason-light` is a lightweight geospatial extension for MASON. It provides support for vector and raster geospatial data, with additional utilities for graph-based spatial modelling, routing, and GIS-oriented agent-based simulations.

The original GeoMason documentation is available from the [GeoMason website](https://cs.gmu.edu/~eclab/projects/mason/extensions/geomason/).

## Versioning note

Versions up to `1.16` followed the original project numbering scheme.

Starting from `2.0.0`, GeoMason-light adopts a clearer versioning scheme. Version `2.0.0` marks the start of the new release line and includes the Java 11 baseline, Maven Central packaging cleanup, package-level Javadocs, and the change of the MASON dependency to `provided`.

## Why GeoMason-light?

### Extended geospatial functionality

`GeoMason-light` builds on the original GeoMason framework and adds a more flexible structure for geospatial and network-based modelling.

Main additions include:

* a `graph` package for advanced manipulation of graph components;
* support for `SubGraph` objects;
* a `VectorLayer` class with additional GIS-style geometric and spatial operations;
* support for loading GeoPackages through [GeoPackage Java](https://github.com/ngageoint/geopackage-java);
* routing utilities for graph-based agent movement and street-network simulations.

These features are intended to support flexible urban mobility, traffic, accessibility, and other spatial agent-based models.

### Maven-based and easier to reuse

The repository is structured as a Maven project, making it easier to build, reuse, document, and integrate into downstream projects.

The library provides:

* Maven-compatible project structure;
* source JAR generation;
* Javadoc JAR generation;
* package-level API documentation;
* CI support for build and Javadoc validation.

### Lightweight distribution

Unlike the original GeoMason distribution, `GeoMason-light` does not include the full set of heavy example models and datasets. The focus is on the reusable core library.

The codebase also uses the modern Java Topology Suite package namespace, `org.locationtech.jts`.

## Requirements

GeoMason-light `2.0.0` requires Java 11 or newer.

## Main dependencies

Core dependencies include:

* `mason`
* `javatuples`
* `jts-core`
* `geopackage-java`

## MASON dependency

`GeoMason-light` depends on MASON.

MASON is not available from Maven Central under the `sim:mason:21` coordinates. For this reason, the dependency is declared as `provided` in the `GeoMason-light` POM.

This means that Maven will not try to pull MASON transitively when users add `GeoMason-light` to their project. Users must provide MASON manually.

If you are building `GeoMason-light` from source, or if your downstream Maven project needs to compile against MASON classes, install `mason-21.jar` into your local Maven repository:

```bash
mvn install:install-file \
	-Dfile=lib/mason-21.jar \
	-DgroupId=sim \
	-DartifactId=mason \
	-Dversion=21 \
	-Dpackaging=jar
```

On Windows PowerShell:

```powershell
mvn install:install-file "-Dfile=.\lib\mason-21.jar" "-DgroupId=sim" "-DartifactId=mason" "-Dversion=21" "-Dpackaging=jar"
```

After that, downstream Maven projects can explicitly declare both dependencies:

```xml
<dependency>
	<groupId>uk.ac.liv.gdsl</groupId>
	<artifactId>GeoMason-light</artifactId>
	<version>2.0.0</version>
</dependency>

<dependency>
	<groupId>sim</groupId>
	<artifactId>mason</artifactId>
	<version>21</version>
</dependency>
```

## Installation

### Option 1: Use GeoMason-light as a Maven dependency

Add the following dependency to your project `pom.xml`:

```xml
<dependency>
	<groupId>uk.ac.liv.gdsl</groupId>
	<artifactId>GeoMason-light</artifactId>
	<version>2.0.0</version>
</dependency>
```

If your project directly uses MASON classes, install `mason-21.jar` locally as described above and add the explicit `sim:mason:21` dependency to your own `pom.xml`.

Then run:

```bash
mvn clean install
```

### Option 2: Use GeoMason-light in Eclipse

1. Make sure Maven support is installed and enabled in Eclipse.
2. Make sure the project uses Java 11 or newer.
3. Create or open a Maven project.
4. Add the `GeoMason-light` dependency to your project `pom.xml`:

```xml
<dependency>
	<groupId>uk.ac.liv.gdsl</groupId>
	<artifactId>GeoMason-light</artifactId>
	<version>2.0.0</version>
</dependency>
```

5. If your project needs MASON classes, install `mason-21.jar` locally and add the explicit `sim:mason:21` dependency.
6. Update the Maven project:

```text
Right click project → Maven → Update Project
```

You can then import and use GeoMason-light classes in your Java code.

## Building GeoMason-light from source

Clone the repository:

```bash
git clone https://github.com/g-filomena/GeoMason-light.git
cd GeoMason-light
```

Install MASON into your local Maven repository:

```bash
mvn install:install-file \
	-Dfile=lib/mason-21.jar \
	-DgroupId=sim \
	-DartifactId=mason \
	-Dversion=21 \
	-Dpackaging=jar
```

On Windows PowerShell:

```powershell
mvn install:install-file "-Dfile=.\lib\mason-21.jar" "-DgroupId=sim" "-DartifactId=mason" "-Dversion=21" "-Dpackaging=jar"
```

Then build the project:

```bash
mvn clean install
```

To generate the Javadoc locally:

```bash
mvn -DskipTests javadoc:javadoc
```

The generated documentation will be available under:

```text
target/site/apidocs/
```

## API documentation

The project includes package-level Javadocs for the main API areas:

* `sim.field.geo`
* `sim.graph`
* `sim.io.geo`
* `sim.portrayal.geo`
* `sim.routing`
* `sim.util.geo`

These packages cover vector/raster layers, geospatial import/export, graph structures, routing, portrayal, and geometry utilities.

## License

GeoMason-light is released under the GNU General Public License v3.0.

See the [GPLv3 license](https://www.gnu.org/licenses/gpl-3.0.en.html) for details.
