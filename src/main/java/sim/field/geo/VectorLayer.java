/*
 * Copyright (c) 2023 Gabriele Filomena
 * University of Liverpool, UK
 * 
 * This program is free software: it can redistributed and/or modified
 * under the terms of the GNU General Public License 3.0 as published by
 * the Free Software Foundation.
 *
 * See the file "LICENSE" for more information
 */
package sim.field.geo;

import java.net.URL;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import org.locationtech.jts.algorithm.ConvexHull;
import org.locationtech.jts.geom.Coordinate;
import org.locationtech.jts.geom.Envelope;
import org.locationtech.jts.geom.Geometry;
import org.locationtech.jts.geom.GeometryFactory;
import org.locationtech.jts.index.quadtree.Quadtree;

import sim.io.geo.ShapeFileImporter;
import sim.util.geo.MasonGeometry;

/**
 * VectorLayer is an extension of GeomVectorField that offers advanced geometric, selection,
 * and filtering capabilities for managing and working with spatial data.
 * It enhances the functionality provided by GeomVectorField.
 */
public class VectorLayer extends GeomVectorField {

	private static final long serialVersionUID = 1L;
	public ArrayList<MasonGeometry> geometriesList = new ArrayList<>();
	private Quadtree layerSpatialIndex = new Quadtree();
	private final GeometryFactory layerGeomFactory = new GeometryFactory();

	/**
	 * Creates an empty VectorLayer with no initial geometries.
	 */
	public VectorLayer() {
	    super();
	}

	/**
	 * Creates a VectorLayer from a list of MasonGeometry objects contained in a Bag.
	 *
	 * @param geometries A Bag containing MasonGeometry objects to populate the VectorLayer with.
	 */
	public VectorLayer(ArrayList<MasonGeometry> geometries) {
		super();

		for (final MasonGeometry masonGeometry : geometries) {
			super.addGeometry(masonGeometry);
			final Envelope envelope = masonGeometry.getGeometry().getEnvelopeInternal();
			layerSpatialIndex.insert(envelope, masonGeometry);
		}
		generateGeometriesList();
	}

	/**
	 * Adds a MasonGeometry to the VectorLayer and updates the layer's spatial index.
	 * This method adds the specified MasonGeometry to the VectorLayer and inserts its envelope
	 * into the layer's spatial index for efficient spatial queries.
	 *
	 * @param masonGeometry The MasonGeometry to be added to the VectorLayer.
	 */
	@Override
	public void addGeometry(MasonGeometry masonGeometry) {
		super.addGeometry(masonGeometry);
		final Envelope envelope = masonGeometry.getGeometry().getEnvelopeInternal();
		layerSpatialIndex.insert(envelope, masonGeometry);
	}

	/**
	 * Retrieves all the geometries in the VectorLayer that intersect with a given Geometry.
	 *
	 * @param inputGeometry The reference geometry defining the intersection criterion.
	 * @return An ArrayList containing the geometries that intersect with the specified input geometry.
	 */
	public final ArrayList<MasonGeometry> intersectingFeatures(Geometry inputGeometry) {

		ArrayList<MasonGeometry> intersectingFeatures = new ArrayList<>();
		final Envelope envelope = inputGeometry.getEnvelopeInternal();
		envelope.expandBy(java.lang.Math.max(envelope.getHeight(), envelope.getWidth()) * 0.01);
		final List<?> geometriesList = layerSpatialIndex.query(envelope);

		for (final Object geometry : geometriesList) {
			final MasonGeometry masonGeometry = (MasonGeometry) geometry;
			if (inputGeometry.intersects(masonGeometry.geometry))
				intersectingFeatures.add(masonGeometry);
		}
		return intersectingFeatures;
	}

	/**
	 * Retrieves all the geometries in the VectorLayer that are fully contained within
	 * the boundaries of a given Geometry.
	 *
	 * @param inputGeometry The reference geometry defining the containment boundary.
	 * @return An ArrayList containing the geometries that are entirely contained within
	 *         the specified input geometry.
	 */
	public final ArrayList<MasonGeometry> containedFeatures(Geometry inputGeometry) {
		final ArrayList<MasonGeometry> containedFeatures = new ArrayList<>();
		final Envelope envelope = inputGeometry.getEnvelopeInternal();
		envelope.expandBy(java.lang.Math.max(envelope.getHeight(), envelope.getWidth()) * 0.01);
		final List<?> geometriesList = layerSpatialIndex.query(envelope);

		for (final Object geometry : geometriesList) {
			final MasonGeometry masonGeometry = (MasonGeometry) geometry;
			if (inputGeometry.contains(masonGeometry.geometry))
				containedFeatures.add(masonGeometry);
		}
		return containedFeatures;
	}

	/**
	 * Retrieves all the geometries in the VectorLayer that are within a specified range of distances
	 * from a given Geometry.
	 *
	 * @param inputGeometry The reference geometry for distance calculation.
	 * @param lowerLimit    The minimum distance from the input geometry (inclusive).
	 * @param upperLimit    The maximum distance from the input geometry (inclusive).
	 * @return An ArrayList containing the geometries that fall within the specified distance range
	 *         from the input geometry.
	 */
	public ArrayList<MasonGeometry> featuresBetweenLimits(Geometry inputGeometry, double lowerLimit,
			double upperLimit) {

		final ArrayList<MasonGeometry> geometries = new ArrayList<>();
		final Envelope envelope = inputGeometry.getEnvelopeInternal();
		envelope.expandBy(upperLimit);
		final List<?> geometriesList = layerSpatialIndex.query(envelope);

		for (final Object geometry : geometriesList) {
			final MasonGeometry masonGeometry = (MasonGeometry) geometry;
			if (inputGeometry.distance(masonGeometry.geometry) >= lowerLimit
					& inputGeometry.distance(masonGeometry.getGeometry()) <= upperLimit)
				geometries.add(masonGeometry);
			else
				continue;
		}
		return geometries;
	}

	/**
	 * Returns all the geometries in the VectorLayer that are contained within a
	 * certain radius from a given Geometry.
	 *
	 * @param inputGeometry The reference geometry for containment evaluation.
	 * @param radius        The maximum distance from the input geometry within which other geometries are considered.
	 * @return An ArrayList containing the geometries that are within the specified radius of the input geometry.
	 */
	public ArrayList<MasonGeometry> featuresWithinDistance(Geometry inputGeometry, double radius) {

		final ArrayList<MasonGeometry> geometries = new ArrayList<>();
		final Envelope envelope = inputGeometry.getEnvelopeInternal();
		envelope.expandBy(radius);
		final List<?> geometriesList = layerSpatialIndex.query(envelope);

		for (final Object geometry : geometriesList) {
			final MasonGeometry masonGeometry = (MasonGeometry) geometry;
			if (inputGeometry.isWithinDistance(masonGeometry.geometry, radius))
				geometries.add(masonGeometry);
		}
		return geometries;
	}

	/**
	 * Filters a list of MasonGeometry objects based on a specified integer attribute's value.
	 * This method iterates through a list of MasonGeometry objects and selects geometries
	 * that either match or do not match the provided integer attribute value based on the 'equal' parameter.
	 * Matching geometries are included in the resulting ArrayList.
	 *
	 * @param attributeName The name of the integer attribute to filter by. This should be a valid attribute
	 *                      name present in the MasonGeometry objects.
	 * @param attributeValue The integer value to compare against for filtering. Geometries with the specified
	 *                      integer attribute equal to this value will be included in the result.
	 * @param equal If true, filter geometries where the specified integer attribute equals the given value.
	 *              If false, filter geometries where the specified integer attribute does not equal the given value.
	 * @return An ArrayList containing the filtered MasonGeometry objects. The list may be empty if
	 *         no matching geometries are found.
	 */
	public ArrayList<MasonGeometry> filterFeatures(String attributeName, int attributeValue, boolean equal) {

		final ArrayList<MasonGeometry> geometries = new ArrayList<>();
		for (final MasonGeometry masonGeometry : geometriesList) {
			final Integer attribute = masonGeometry.getIntegerAttribute(attributeName);
			if (!equal && !attribute.equals(attributeValue))
				geometries.add(masonGeometry);
			else if (attribute.equals(attributeValue))
				geometries.add(masonGeometry);
		}
		return geometries;
	}

	/**
	 * Filters a list of MasonGeometry objects based on a specified attribute's value.
	 * This method iterates through a list of MasonGeometry objects and selects geometries
	 * that either match or do not match the provided attribute value based on the 'equal' parameter.
	 * Matching geometries are included in the resulting ArrayList.
	 *
	 * @param attributeName The name of the attribute to filter by. This should be a valid attribute
	 *                      name present in the MasonGeometry objects.
	 * @param attributeValue The value to compare against for filtering. Geometries with the specified
	 *                      attribute equal to this value will be included in the result.
	 * @param equal If true, filter geometries where the specified attribute equals the given value.
	 *              If false, filter geometries where the specified attribute does not equal the given value.
	 * @return An ArrayList containing the filtered MasonGeometry objects. The list may be empty if
	 *         no matching geometries are found
	 */
	public ArrayList<MasonGeometry> filterFeatures(String attributeName, String attributeValue, boolean equal) {

		final ArrayList<MasonGeometry> geometries = new ArrayList<>();
		for (final MasonGeometry masonGeometry : geometriesList) {
			final String attribute = masonGeometry.getStringAttribute(attributeName);
			if (!equal && !attribute.equals(attributeValue))
				geometries.add(masonGeometry);
			else if (attribute.equals(attributeValue))
				geometries.add(masonGeometry);
		}
		return geometries;
	}

	/**
	 * Filters and returns a list of MasonGeometry objects based on a specified string attribute and values.
	 * This function filters the features in this VectorLayer based on the provided string attribute and a list of
	 * string values. If 'equal' is true, it keeps features with attribute values that match those in the list;
	 * if 'equal' is false, it keeps features with attribute values that do not match the list. The selected features
	 * are added to a new list and returned.
	 *
	 * @param attributeName The name of the string attribute to filter on.
	 * @param listValues    A list of string values used for filtering.
	 * @param equal         If true, keeps features with attribute values matching the list; if false, keeps others.
	 * @return              A list of MasonGeometry objects after filtering.
	 */
	public ArrayList<MasonGeometry> filterFeatures(String attributeName, List<String> listValues, boolean equal) {

		final ArrayList<MasonGeometry> geometries = new ArrayList<>();
		for (final MasonGeometry masonGeometry : geometriesList) {
			final String attribute = masonGeometry.getStringAttribute(attributeName);
			if (!equal && !listValues.contains(attribute))
				geometries.add(masonGeometry);
			else if (listValues.contains(attribute))
				geometries.add(masonGeometry);
		}
		return null;
	}

	/**
	 * Retrieves a list of integer values from a specified attribute column in this VectorLayer.
	 * This function iterates through the geometries in the VectorLayer and extracts integer values
	 * from the specified attribute column for each geometry. It returns a list of these integer values.
	 *
	 * @param attributeName The name of the attribute column to retrieve values from.
	 * @return              A list of integer values from the specified attribute column.
	 */
	public List<Integer> getIntColumn(String attributeName) {

		final List<Integer> values = new ArrayList<>();
		for (final MasonGeometry masonGeometry : geometriesList) {
			final Integer attribute = masonGeometry.getIntegerAttribute(attributeName);
			values.add(attribute);
		}
		return values;
	}

	/**
	 * Selects and creates a new VectorLayer containing features based on a specified integer attribute and values.
	 * This function filters the features in this VectorLayer based on the provided integer attribute and a list of
	 * integer values. If 'equal' is true, it keeps features with attribute values that match those in the list;
	 * if 'equal' is false, it keeps features with attribute values that do not match the list. The selected features
	 * are added to a new VectorLayer, which is returned.
	 *
	 * @param attributeName The name of the integer attribute to filter on.
	 * @param listValues    A list of integer values used for filtering.
	 * @param equal         If true, keeps features with attribute values matching the list; if false, keeps others.
	 * @return              A new VectorLayer containing the selected features.
	 */
	public VectorLayer selectFeatures(String attributeName, List<Integer> listValues, boolean equal) {

		ArrayList<MasonGeometry> geometries = new ArrayList<>();

		for (final MasonGeometry masonGeometry : geometriesList) {
			final Integer attribute = masonGeometry.getIntegerAttribute(attributeName);
			if (!equal && !listValues.contains(attribute))
				geometries.add(masonGeometry);
			else if (equal && listValues.contains(attribute))
				geometries.add(masonGeometry);
		}
		final VectorLayer newLayer = new VectorLayer(geometries);
		return newLayer;
	}

	/**
	 * Computes the convex hull of the geometries within this VectorLayer.
	 * This function calculates the convex hull, which is the smallest convex polygon that encloses
	 * all the geometries in this VectorLayer. It returns the convex hull as a Geometry object.
	 *
	 * @return The convex hull of the geometries within this VectorLayer.
	 */
	public Geometry layerConvexHull() {
		final ArrayList<Coordinate> pts = new ArrayList<>();

		for (final MasonGeometry masonGeometry : geometriesList) {
			final Geometry geometry = masonGeometry.geometry;
			final Coordinate coordinateGeometry[] = geometry.getCoordinates();
			pts.addAll(Arrays.asList(coordinateGeometry));
		}

		final Coordinate[] coords = pts.toArray(new Coordinate[pts.size()]);
		final ConvexHull convexHull = new ConvexHull(coords, layerGeomFactory);
		return convexHull.getConvexHull();
	}

	/**
	 * Retrieves a list of integer IDs associated with the geometries in this VectorLayer.
	 * This function iterates through the geometries in the VectorLayer and extracts the integer IDs
	 * stored as user data in each MasonGeometry. It returns a list of these IDs.
	 *
	 * @return A list of integer IDs associated with the geometries in this VectorLayer.
	 */
	public void setID(String attributeName) {

		for (final MasonGeometry masonGeometry : geometriesList)
			masonGeometry.setUserData(masonGeometry.getIntegerAttribute(attributeName));
	}

	/**
	 * It gets the ID of this laye
	 *
	 * @return
	 */
	public ArrayList<Integer> getIDs() {

		ArrayList<Integer> IDs = new ArrayList<>();
		for (final MasonGeometry masonGeometry : geometriesList)
			IDs.add((Integer) masonGeometry.getUserData());
		return IDs;
	}

	/**
	 * Computes the intersection of geometries between this VectorLayer and another VectorLayer.
	 * This function takes another VectorLayer as input and checks for intersections between its geometries
	 * and the geometries in this VectorLayer. It returns a list of intersecting geometries if 'inclusive' is true,
	 * or a list of non-intersecting geometries if 'inclusive' is false.
	 *
	 * @param otherLayer  The VectorLayer to intersect with this VectorLayer.
	 * @param inclusive   If true, returns intersecting geometries; if false, returns non-intersecting geometries.
	 * @return            A list of intersecting or non-intersecting MasonGeometry objects.
	 */
	public ArrayList<MasonGeometry> intersection(VectorLayer otherLayer, boolean inclusive) {

		ArrayList<MasonGeometry> intersectingGeometries = new ArrayList<>();
		for (final MasonGeometry masonGeometry : otherLayer.geometriesList)
			intersectingGeometries.addAll(intersectingFeatures(masonGeometry.geometry));
		if (inclusive)
			return intersectingGeometries;
		else {
			ArrayList<MasonGeometry> notIntersecting = otherLayer.geometriesList;
			notIntersecting.removeAll(intersectingGeometries);
			return notIntersecting;
		}
	}

	/**
	 * Verifies if any of the features of this VectorLayer intersect a given
	 * geometry.
	 *
	 * @param inputGeometry the input geometry on which the intersection is
	 *                      verified;
	 */
	public boolean intersects(Geometry inputGeometry) {

		final Envelope envelope = inputGeometry.getEnvelopeInternal();
		envelope.expandBy(java.lang.Math.max(envelope.getHeight(), envelope.getWidth()) * 0.01);
		final List<?> geometriesList = layerSpatialIndex.query(envelope);

		for (final Object geometry : geometriesList) {
			final MasonGeometry masonGeometry = (MasonGeometry) geometry;
			if (inputGeometry.intersects(masonGeometry.geometry))
				return true;
		}
		return false;
	}

	/**
	 * Generates a list of MasonGeometry objects from the geometries in the VectorLayer.
	 * This function clears the existing geometries list, iterates through the geometries in the VectorLayer,
	 * adds each geometry to the list, and inserts their envelopes into the spatial index of the layer.
	 * It is used to refresh the geometries list after changes in the VectorLayer.
	 */
	public void generateGeometriesList() {

		geometriesList.clear();
		layerSpatialIndex = new Quadtree();
		for (final Object geometry : getGeometries()) {
			final MasonGeometry masonGeometry = (MasonGeometry) geometry;
			geometriesList.add(masonGeometry);
			final Envelope envelope = masonGeometry.getGeometry().getEnvelopeInternal();
			layerSpatialIndex.insert(envelope, masonGeometry);
		}
	}
	
	/**
	 * Reads a Shapefile from the specified input directory and populates the provided VectorLayer with geometries.
	 * The Shapefile consists of two files: a .shp file containing geometry information and a .dbf file containing attribute data.
	 * The function reads both files and adds the geometries to the VectorLayer. After reading, it generates a list of geometries
	 * for the VectorLayer.
	 *
	 * @param inputDirectory The directory path where the Shapefile (.shp and .dbf) is located.
	 * @param vectorLayer    The VectorLayer to populate with the geometries from the Shapefile.
	 * @throws Exception     If there is an error during Shapefile reading or geometry population.
	 */
	public static void readShapefile(String inputDirectory, VectorLayer vectorLayer) throws Exception {

		String baseURL = "file:";
		final URL shp = new URL(baseURL + inputDirectory + ".shp");
		final URL dbf = new URL(baseURL + inputDirectory + ".dbf");
		ShapeFileImporter.read(shp, dbf, vectorLayer);
		vectorLayer.generateGeometriesList();
	}
}