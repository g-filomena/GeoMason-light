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

import java.awt.geom.AffineTransform;
import java.awt.geom.Point2D;
import java.net.URL;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.concurrent.ConcurrentLinkedQueue;

import org.locationtech.jts.algorithm.ConvexHull;
import org.locationtech.jts.geom.Coordinate;
import org.locationtech.jts.geom.CoordinateSequenceFilter;
import org.locationtech.jts.geom.Envelope;
import org.locationtech.jts.geom.Geometry;
import org.locationtech.jts.geom.GeometryFactory;
import org.locationtech.jts.geom.Point;
import org.locationtech.jts.geom.Polygon;
import org.locationtech.jts.geom.prep.PreparedGeometryFactory;
import org.locationtech.jts.geom.prep.PreparedPolygon;
import org.locationtech.jts.index.quadtree.Quadtree;

import sim.engine.SimState;
import sim.engine.Steppable;
import sim.io.geo.GeoPackageImporter;
import sim.io.geo.ShapeFileImporter;
import sim.portrayal.DrawInfo2D;
import sim.util.geo.GeometryUtilities;
import sim.util.geo.MasonGeometry;

/**
 * `VectorLayer` is an extension of the `Layer` class that offers advanced
 * geometric, selection, and filtering capabilities for managing and working
 * with spatial data.
 */
public class VectorLayer extends Layer {

	private static final long serialVersionUID = 1L;
	private List<MasonGeometry> geometriesList = new ArrayList<>();
	private final GeometryFactory layerGeomFactory = new GeometryFactory();

	/**
	 * A spatial index of all the geometries in the field.
	 */
	private Quadtree spatialIndex = new Quadtree();

	/**
	 * The convex hull of all the geometries in this field.
	 */
	private PreparedPolygon convexHull = null;

	/**
	 * Defines the outer shell of all the geometries within this field.
	 */
	private PreparedPolygon union;

	/**
	 * Helper factory for computing the union or convex hull.
	 */
	private GeometryFactory geomFactory = new GeometryFactory();

	// Constructors

	/**
	 * Creates an empty VectorLayer.
	 */
	public VectorLayer() {
		super();
	}

	/**
	 * Creates an empty VectorLayer with no initial geometries.
	 * 
	 * @param width  The field width in display units for managing scale changes.
	 * @param height The field height in display units for managing scale changes.
	 */
	public VectorLayer(int width, int height) {
		super(width, height);
	}

	/**
	 * Creates a VectorLayer from a List of MasonGeometry objects contained in an
	 * List of MasonGeometries.
	 *
	 * @param masonGeometries A List containing MasonGeometry objects to populate
	 *                        the VectorLayer with.
	 */
	public VectorLayer(List<MasonGeometry> masonGeometries) {
		super();
		MBR = new Envelope();
		for (final MasonGeometry masonGeometry : masonGeometries)
			addGeometry(masonGeometry);
	}

	// Basic Setters/Getters

	/**
	 * Adds a MasonGeometry to the VectorLayer and updates the layer's spatial
	 * index. This method adds the specified MasonGeometry to the VectorLayer and
	 * inserts its envelope into the layer's spatial index for efficient spatial
	 * queries.
	 *
	 * @param masonGeometry The MasonGeometry to be added to the VectorLayer.
	 */
	public void addGeometry(MasonGeometry masonGeometry) {

		final Envelope envelope = masonGeometry.getGeometry().getEnvelopeInternal();
		MBR.expandToInclude(envelope);
		spatialIndex.insert(envelope, masonGeometry);
		geometriesList.add(masonGeometry);
	}

	/**
	 * Searches the VectorLayer for the first geometry with attribute equals to name
	 * that has the given value.
	 * 
	 * @param name  of attribute.
	 * @param value of attribute.
	 * @return MasonGeometry with specified attribute otherwise null.
	 */
	public MasonGeometry getGeometry(String name, Object value) {

		for (MasonGeometry masonGeometry : geometriesList) {
			if (masonGeometry.hasAttribute(name) && masonGeometry.getAttribute(name).equals(value)) {
				return masonGeometry;
			}
		}
		return null;
	}

	/**
	 * Removes the given geometry.
	 * 
	 * @param masonGeometry The MasonGeometry to be removed from the VectorLayer.
	 */
	public void removeGeometry(final MasonGeometry masonGeometry) {
		geometriesList.remove(masonGeometry);
	}

	/**
	 * Returns the List of MasonGeometry features.
	 *
	 * @return The List of MasonGeometry features.
	 */
	public List<MasonGeometry> getGeometries() {
		return new ArrayList<>(geometriesList);
	}

	/**
	 * Sets a unique identifier for each MasonGeometry object in the geometries
	 * list.
	 *
	 * This method assigns a unique identifier to each MasonGeometry object in the
	 * geometries list based on the specified attribute name. The attribute value
	 * for each MasonGeometry is used as its unique identifier.
	 *
	 * @param attributeName The name of the attribute that contains unique
	 *                      identifiers.
	 */
	public void setID(String attributeName) {
		for (final MasonGeometry masonGeometry : geometriesList)
			masonGeometry.setUserData(masonGeometry.getIntegerAttribute(attributeName));
	}

	/**
	 * Retrieves a List of IDs from the MasonGeometry objects in the geometries
	 * list. Each ID is obtained from the user data of the MasonGeometry object,
	 * which is assumed to be stored as an Integer. This method iterates through all
	 * MasonGeometry objects in the geometriesList and collects their IDs into a new
	 * list.
	 *
	 * @return A List of Integer IDs, one for each MasonGeometry object in the
	 *         geometries list.
	 */
	public List<Integer> getIDs() {
		List<Integer> IDs = new ArrayList<>();
		for (final MasonGeometry masonGeometry : geometriesList)
			IDs.add((Integer) masonGeometry.getUserData()); // this needed to be set
		return IDs;
	}

	/**
	 * Retrieves a List of integer values from a specified attribute column in this
	 * VectorLayer. This function iterates through the geometries in the VectorLayer
	 * and extracts integer values from the specified attribute column for each
	 * geometry. It returns a List of these integer values.
	 *
	 * @param attributeName The name of the attribute column to retrieve values
	 *                      from.
	 * @return A List of integer values from the specified attribute column.
	 */
	public List<Integer> getIntColumn(String attributeName) {
		final List<Integer> values = new ArrayList<>();
		for (final MasonGeometry masonGeometry : geometriesList) {
			final Integer attribute = masonGeometry.getIntegerAttribute(attributeName);
			values.add(attribute);
		}
		return values;
	}

	// Geometry Manipulation

	/**
	 * Moves the centroid of the given geometry to the provided point. The spatial
	 * index is not notified of the geometry changes. It is strongly recommended
	 * that updateSpatialIndex() be invoked after all geometry position changes.
	 *
	 * @param masonGeometry The geometry to move.
	 * @param coordsFilter  The coordinate sequence filter.
	 */
	public void setGeometryLocation(MasonGeometry masonGeometry, CoordinateSequenceFilter coordsFilter) {
		MasonGeometry otherMasonGeometry = findGeometry(masonGeometry);
		if (otherMasonGeometry != null) {
			otherMasonGeometry.geometry.apply(coordsFilter);
			otherMasonGeometry.geometry.geometryChanged();
		}
		updateSpatialIndex();
	}

	/**
	 * Sets the geometry location of the given MasonGeometry to a new location and
	 * updates the spatial index.
	 *
	 * @param masonGeometry the MasonGeometry object to update
	 * @param newLocation   the new location to set for the MasonGeometry
	 */
	public void setGeometryLocation(MasonGeometry masonGeometry, Point newLocation) {
		MasonGeometry otherMasonGeometry = findGeometry(masonGeometry);
		otherMasonGeometry.geometry = newLocation;
		updateSpatialIndex();
	}

	/**
	 * Retrieves the location of the geometry as a Point object.
	 *
	 * @param masonGeometry The geometry for which the location is to be retrieved.
	 * @return The centroid of the geometry as a Point object.
	 */
	public Point getGeometryLocation(MasonGeometry masonGeometry) {
		MasonGeometry otherMasonGeometry = findGeometry(masonGeometry);
		if (otherMasonGeometry.equals(masonGeometry)) {
			return otherMasonGeometry.getGeometry().getCentroid();
		}
		return null;
	}

	/**
	 * Finds a MasonGeometry within the spatial index.
	 *
	 * @param masonGeometry The geometry to be found.
	 * @return The geometry if found, otherwise returns the original input geometry.
	 */
	public synchronized MasonGeometry findGeometry(MasonGeometry masonGeometry) {

		List<?> geometriesList = spatialIndex.query(masonGeometry.getGeometry().getEnvelopeInternal());

		for (final Object geometry : geometriesList) {
			final MasonGeometry otherMasonGeometry = (MasonGeometry) geometry;
			if (otherMasonGeometry.equals(masonGeometry)) {
				return otherMasonGeometry;
			}
		}
		return masonGeometry;
	}

	// Spatial Index

	/**
	 * Updates the spatial index with the current geometries in the geometriesList.
	 */
	public void updateSpatialIndex() {
		spatialIndex = new Quadtree();
		for (MasonGeometry masonGeometry : geometriesList)
			spatialIndex.insert(masonGeometry.getGeometry().getEnvelopeInternal(), masonGeometry);
	}

	/**
	 * Schedules an updater for the spatial index.
	 *
	 * @return A Steppable object that updates the spatial index
	 */
	public Steppable scheduleSpatialIndexUpdater() {
		return new Steppable() {
			private static final long serialVersionUID = 1808772743124538274L;

			@Override
			public void step(SimState state) {
				updateSpatialIndex();
			}
		};
	}

	/**
	 * Removes all the geometries and resets the MBR.
	 */
	@Override
	public void clear() {
		super.clear();
		spatialIndex = new Quadtree();
		geometriesList.clear();
	}

	// Geometric Computations

	/**
	 * Computes the convex hull of the geometries within this VectorLayer. This
	 * function calculates the convex hull, which is the smallest convex polygon
	 * that encloses all the geometries in this VectorLayer.
	 */
	private void computeConvexHull() {
		final List<Coordinate> points = new ArrayList<>();

		if (geometriesList.isEmpty()) {
			return;
		}

		for (final MasonGeometry masonGeometry : geometriesList) {
			final Geometry geometry = masonGeometry.getGeometry();
			final Coordinate coordinateGeometry[] = geometry.getCoordinates();
			points.addAll(Arrays.asList(coordinateGeometry));
		}

		final Coordinate[] coordinates = points.toArray(new Coordinate[points.size()]);
		final ConvexHull convexHull = new ConvexHull(coordinates, layerGeomFactory);
		this.convexHull = new PreparedPolygon((Polygon) convexHull.getConvexHull());
	}

	/**
	 * Computes and retrieves the convex hull of the geometries within this
	 * VectorLayer. The convex hull is the smallest convex polygon that encloses all
	 * the geometries in this VectorLayer.
	 *
	 * @return The convex hull geometry of the geometries in this VectorLayer
	 */
	public Geometry getConvexHull() {
		if (convexHull == null)
			computeConvexHull();
		return convexHull.getGeometry();
	}

	/**
	 * Compute the union of the VectorLayer geometries. The resulting Geometry is
	 * the outside points of the layer geometries.
	 */
	private void computeUnion() {
		Geometry polygon = new Polygon(null, null, geomFactory);

		if (geometriesList.isEmpty())
			return;

		for (MasonGeometry masonGeometry : geometriesList) {
			Geometry geometry = masonGeometry.getGeometry();
			polygon = polygon.union(geometry);
		}

		polygon = polygon.union();
		union = new PreparedPolygon((Polygon) polygon);
	}

	/**
	 * Computes and retrieves the union of the geometries within this VectorLayer.
	 * The resulting Geometry represents the outside points of the field's
	 * geometries.
	 *
	 * @return The union geometry of the geometries in this VectorLayer.
	 */
	public Geometry getUnion() {
		if (union == null) {
			computeUnion();
		}
		return union.getGeometry();
	}

	// Spatial Queries

	/**
	 * Queries the VectorLayer based on the provided Envelope and returns a List of
	 * MasonGeometry objects that intersect with the envelope.
	 *
	 * @param envelope The envelope used for the query.
	 * @return A List of MasonGeometry objects that intersect with the provided
	 *         envelope.
	 */
	public synchronized List<MasonGeometry> queryField(Envelope envelope) {
		List<?> geometriesList = spatialIndex.query(envelope);
		List<MasonGeometry> geometries = new ArrayList<MasonGeometry>(geometriesList.size());

		// However, the JTS QuadTree query is a little sloppy, which means it
		// may return objects that are still outside the range. We need to do
		// a second pass to trim out the objects that are further than distance.
		for (final Object geometry : geometriesList) {
			MasonGeometry masonGeometry = (MasonGeometry) geometry;
			if (envelope.intersects(masonGeometry.getGeometry().getEnvelopeInternal())) {
				geometries.add(masonGeometry);
			}
		}
		return geometries;
	}

	/**
	 * Retrieves all the geometries in the VectorLayer that are within a specified
	 * range of distances from a given Geometry.
	 *
	 * @param inputGeometry The reference geometry for distance calculation.
	 * @param lowerLimit    The minimum distance from the input geometry
	 *                      (inclusive).
	 * @param upperLimit    The maximum distance from the input geometry
	 *                      (inclusive).
	 * @return A List containing the geometries that fall within the specified
	 *         distance range from the input geometry.
	 */
	public List<MasonGeometry> featuresBetweenLimits(Geometry inputGeometry, double lowerLimit, double upperLimit) {
		final List<MasonGeometry> featuresBetween = new ArrayList<>();
		final Envelope envelope = inputGeometry.getEnvelopeInternal();
		envelope.expandBy(upperLimit);
		final List<?> geometriesList = spatialIndex.query(envelope);

		for (final Object geometry : geometriesList) {
			final MasonGeometry otherMasonGeometry = (MasonGeometry) geometry;
			if (inputGeometry.distance(otherMasonGeometry.getGeometry()) >= lowerLimit
					& inputGeometry.distance(otherMasonGeometry.getGeometry()) <= upperLimit)
				featuresBetween.add(otherMasonGeometry);
			else
				continue;
		}
		return featuresBetween;
	}

	/**
	 * Returns all the geometries in the VectorLayer that are contained within a
	 * certain radius from a given Geometry.
	 *
	 * @param inputGeometry The reference geometry for containment evaluation.
	 * @param radius        The maximum distance from the input geometry within
	 *                      which other geometries are considered.
	 * @return A List containing the geometries that are within the specified radius
	 *         of the input geometry.
	 */
	public List<MasonGeometry> featuresWithinDistance(Geometry inputGeometry, double radius) {
		final List<MasonGeometry> featuresWithin = new ArrayList<>();
		final Envelope envelope = inputGeometry.getEnvelopeInternal();
		envelope.expandBy(radius);
		final List<?> geometriesList = spatialIndex.query(envelope);

		for (final Object geometry : geometriesList) {
			MasonGeometry otherMasonGeometry = (MasonGeometry) geometry;
			if (inputGeometry.isWithinDistance(otherMasonGeometry.getGeometry(), radius))
				featuresWithin.add(otherMasonGeometry);
		}
		return featuresWithin;
	}

	/**
	 * Retrieves the geometries in the VectorLayer that intersect with a given
	 * Geometry.
	 *
	 * @param inputGeometry The reference geometry defining the intersection
	 *                      criterion.
	 * @return A List containing the geometries that intersect with the specified
	 *         input geometry.
	 */
	public final List<MasonGeometry> intersectingFeatures(Geometry inputGeometry) {
		ConcurrentLinkedQueue<MasonGeometry> intersectingFeatures = new ConcurrentLinkedQueue<>();
		final Envelope envelope = inputGeometry.getEnvelopeInternal();
		envelope.expandBy(Math.max(envelope.getHeight(), envelope.getWidth()) * 0.01);
		final List<?> geometriesList = spatialIndex.query(envelope);

		geometriesList.parallelStream().map(geometry -> (MasonGeometry) geometry)
				.filter(otherMasonGeometry -> inputGeometry.intersects(otherMasonGeometry.getGeometry()))
				.forEach(intersectingFeatures::add);

		return new ArrayList<>(intersectingFeatures);
	}

	/**
	 * Computes the intersection of geometries between this VectorLayer and another
	 * VectorLayer. This function takes another VectorLayer as input and checks for
	 * intersections between its geometries and the geometries in this VectorLayer.
	 * It returns a List of intersecting geometries if 'inclusive' is true, or a
	 * list of non-intersecting geometries if 'inclusive' is false.
	 *
	 * @param otherLayer The VectorLayer to intersect with this VectorLayer.
	 * @param inclusive  If true, returns intersecting geometries; if false, returns
	 *                   non-intersecting geometries.
	 * @return A List of intersecting or non-intersecting MasonGeometry features.
	 */
	public List<MasonGeometry> intersection(VectorLayer otherLayer, boolean inclusive) {
		List<MasonGeometry> intersectingGeometries = new ArrayList<>();
		for (final MasonGeometry masonGeometry : otherLayer.geometriesList)
			intersectingGeometries.addAll(intersectingFeatures(masonGeometry.getGeometry()));
		if (inclusive)
			return intersectingGeometries;
		else {
			List<MasonGeometry> notIntersecting = otherLayer.geometriesList;
			notIntersecting.removeAll(intersectingGeometries);
			return notIntersecting;
		}
	}

	// Spatial Relations

	/**
	 * Checks whether the input MasonGeometry is covered by other geometries in the
	 * spatial index of the VectorLayer.
	 *
	 * @param inputMasonGeometry The MasonGeometry object to check for coverage.
	 * @return true if the inputMasonGeometry is covered by other geometries in the
	 *         spatial index, false otherwise.
	 */
	public boolean isCovered(MasonGeometry inputMasonGeometry) {
		Envelope envelope = inputMasonGeometry.getGeometry().getEnvelopeInternal();
		List<?> geometriesList = spatialIndex.query(envelope);
		if (inputMasonGeometry.preparedGeometry == null) {
			inputMasonGeometry.preparedGeometry = PreparedGeometryFactory.prepare(inputMasonGeometry.getGeometry());
		}

		for (final Object geometry : geometriesList) {
			MasonGeometry otherMasonGeometry = (MasonGeometry) geometry;
			if (!inputMasonGeometry.equals(otherMasonGeometry)
					&& otherMasonGeometry.getGeometry().covers(inputMasonGeometry.getGeometry())) {
				return true;
			}
		}
		return false;
	}

	/**
	 * Checks if the given coordinate is inside the convex hull of the VectorLayer.
	 *
	 * @param coord The coordinate to check.
	 * @return true if the coordinate is inside the convex hull of the VectorLayer,
	 *         false otherwise.
	 */
	public boolean isInsideConvexHull(final Coordinate coord) {
		Point point = geomFactory.createPoint(coord);
		if (convexHull == null) {
			computeConvexHull();
		}
		if (convexHull.intersects(point)) {
			return true;
		}
		return false;
	}

	/**
	 * Checks if the given coordinate is inside the union of the VectorLayer.
	 *
	 * @param coordinate The coordinate to check.
	 * @return true if the coordinate is inside the union of the VectorLayer, false
	 *         otherwise.
	 */
	public boolean isInsideUnion(final Coordinate coordinate) {
		Point point = geomFactory.createPoint(coordinate);
		if (union == null) {
			computeUnion();
		}
		if (union.intersects(point)) {
			return true;
		}
		return false;
	}

	/**
	 * Verifies if any of the geometries of this VectorLayer intersect a given
	 * geometry.
	 *
	 * @param inputGeometry The input geometry on which the intersection is
	 *                      verified.
	 */
	public boolean intersects(Geometry inputGeometry) {
		final Envelope envelope = inputGeometry.getEnvelopeInternal();
		envelope.expandBy(Math.max(envelope.getHeight(), envelope.getWidth()) * 0.01);
		final List<?> geometriesList = spatialIndex.query(envelope);

		for (final Object geometry : geometriesList) {
			final MasonGeometry otherMasonGeometry = (MasonGeometry) geometry;
			if (inputGeometry.intersects(otherMasonGeometry.getGeometry())) {
				return true;
			}
		}
		return false;
	}

	// Feature Filtering

	/**
	 * Filters a List of geometries based on a specified String attribute's value.
	 * This method selects geometries that either match or do not match the provided
	 * attribute value based on the 'equal' parameter. Matching geometries are
	 * included in the resulting List.
	 *
	 * @param attributeName  The name of the attribute to filter by. This should be
	 *                       a valid attribute name present in the MasonGeometry
	 *                       objects.
	 * @param attributeValue The value to compare against for filtering. Geometries
	 *                       with the specified attribute equal to this value will
	 *                       be included in the result.
	 * @param equal          If true, filter geometries where the specified
	 *                       attribute equals the given value. If false, filter
	 *                       geometries where the specified attribute does not equal
	 *                       the given value.
	 * @return A List containing the filtered MasonGeometry features. The List may
	 *         be empty if no matching geometries are found.
	 */
	public List<MasonGeometry> filterFeatures(String attributeName, String attributeValue, boolean equal) {
		final List<MasonGeometry> filteredFeatures = new ArrayList<>();
		for (final MasonGeometry masonGeometry : geometriesList) {
			final String attribute = masonGeometry.getStringAttribute(attributeName);
			if (!equal && !attribute.equals(attributeValue))
				filteredFeatures.add(masonGeometry);
			else if (attribute.equals(attributeValue))
				filteredFeatures.add(masonGeometry);
		}
		return filteredFeatures;
	}

	/**
	 * Filters a List of geometries based on a List of String values. This function
	 * filters the geometries in this VectorLayer based on the provided string
	 * attribute and a List of string values. If 'equal' is true, it keeps features
	 * with attribute values that match those in the list; if 'equal' is false, it
	 * keeps features with attribute values that do not match the list. The selected
	 * features are added to a new list and returned.
	 *
	 * @param attributeName The name of the string attribute to filter on.
	 * @param listValues    A List of string values used for filtering.
	 * @param equal         If true, keeps features with attribute values matching
	 *                      the list; if false, keeps others.
	 * @return A List of MasonGeometry features after filtering.
	 */
	public List<MasonGeometry> filterFeatures(String attributeName, List<String> listValues, boolean equal) {
		final List<MasonGeometry> filteredFeatures = new ArrayList<>();
		for (final MasonGeometry masonGeometry : geometriesList) {
			final String attribute = masonGeometry.getStringAttribute(attributeName);
			if (!equal && !listValues.contains(attribute))
				filteredFeatures.add(masonGeometry);
			else if (listValues.contains(attribute))
				filteredFeatures.add(masonGeometry);
		}
		return null;
	}

	/**
	 * Filters a List of geometries in the VectorLayer based on a specified integer
	 * attribute's value. This method selects geometries that either match or do not
	 * match the provided integer attribute value based on the 'equal' parameter.
	 * Matching geometries are included in the resulting List.
	 *
	 * @param attributeName  The name of the integer attribute to filter by. This
	 *                       should be a valid attribute name present in the
	 *                       MasonGeometry features.
	 * @param attributeValue The integer value to compare against for filtering.
	 *                       Geometries with the specified integer attribute equal
	 *                       to this value will be included in the result.
	 * @param equal          If true, filter geometries where the specified integer
	 *                       attribute equals the given value. If false, filter
	 *                       geometries where the specified integer attribute does
	 *                       not equal the given value.
	 * @return A List containing the filtered MasonGeometry features. The List may
	 *         be empty if no matching geometries are found.
	 */
	public List<MasonGeometry> filterFeatures(String attributeName, int attributeValue, boolean equal) {

		final List<MasonGeometry> filteredFeatures = new ArrayList<>();
		for (final MasonGeometry masonGeometry : geometriesList) {
			final Integer attribute = masonGeometry.getIntegerAttribute(attributeName);
			if (!equal && !attribute.equals(attributeValue))
				filteredFeatures.add(masonGeometry);
			else if (attribute.equals(attributeValue))
				filteredFeatures.add(masonGeometry);
		}
		return filteredFeatures;
	}

	/**
	 * Selects and creates a new VectorLayer containing features based on a
	 * specified integer attribute and values. This function filters the features in
	 * this VectorLayer based on the provided integer attribute and a List of
	 * integer values. If 'equal' is true, it keeps geometries with attribute values
	 * that match those in the list; if 'equal' is false, it keeps geometries with
	 * attribute values that do not match the list. The selected geometries are
	 * added to a new VectorLayer, which is returned.
	 *
	 * @param attributeName The name of the integer attribute to filter on.
	 * @param listValues    A List of integer values used for filtering.
	 * @param equal         If true, keeps features with attribute values matching
	 *                      the list; if false, keeps others.
	 * @return A new VectorLayer containing the selected features.
	 */
	public VectorLayer selectFeatures(String attributeName, List<Integer> listValues, boolean equal) {
		List<MasonGeometry> selectedFeatures = new ArrayList<>();
		for (final MasonGeometry masonGeometry : geometriesList) {
			final Integer attribute = masonGeometry.getIntegerAttribute(attributeName);
			if (!equal && !listValues.contains(attribute))
				selectedFeatures.add(masonGeometry);
			else if (equal && listValues.contains(attribute))
				selectedFeatures.add(masonGeometry);
		}
		final VectorLayer newLayer = new VectorLayer(selectedFeatures);
		return newLayer;
	}

	// Feature Relationships

	/**
	 * Identifies the geometries in the VectorLayer that cover the input
	 * MasonGeometry.
	 *
	 * @param inputMasonGeometry The MasonGeometry to find covering features for.
	 * @return A List of MasonGeometry objects in the VectorLayer that cover the
	 *         input MasonGeometry.
	 */
	public final List<MasonGeometry> coveringFeatures(MasonGeometry inputMasonGeometry) {
		List<MasonGeometry> coveringFeatures = new ArrayList<MasonGeometry>();
		Envelope envelope = inputMasonGeometry.getGeometry().getEnvelopeInternal();
		List<?> geometriesList = spatialIndex.query(envelope);

		for (final Object geometry : geometriesList) {
			MasonGeometry otherMasonGeometry = (MasonGeometry) geometry;
			if (!inputMasonGeometry.equals(otherMasonGeometry)
					&& otherMasonGeometry.preparedGeometry.covers(inputMasonGeometry.getGeometry()))
				coveringFeatures.add(otherMasonGeometry);
		}
		return coveringFeatures;
	}

	/**
	 * Identifies the geometries in the VectorLayer that are covered by the input
	 * MasonGeometry.
	 *
	 * @param inputMasonGeometry The MasonGeometry to find covered features for.
	 * @return A List of MasonGeometry objects in the VectorLayer that are covered
	 *         by the input MasonGeometry.
	 */
	public final List<MasonGeometry> coveredFeatures(MasonGeometry inputMasonGeometry) {
		List<MasonGeometry> coveredFeatures = new ArrayList<MasonGeometry>();
		if (inputMasonGeometry.preparedGeometry == null) {
			inputMasonGeometry.preparedGeometry = PreparedGeometryFactory.prepare(inputMasonGeometry.getGeometry());
		}

		for (MasonGeometry otherMasonGeometry : geometriesList) {
			Geometry otherGeometry = otherMasonGeometry.getGeometry();
			if (!inputMasonGeometry.equals(otherMasonGeometry)
					&& inputMasonGeometry.preparedGeometry.covers(otherGeometry))
				coveredFeatures.add(otherMasonGeometry);
		}
		return coveredFeatures;
	}

	/**
	 * Identifies the geometries in the VectorLayer that are covered by the input
	 * MasonGeometry.
	 *
	 * @param inputMasonGeometry The MasonGeometry to find covered features for.
	 * @return A List of MasonGeometry objects in the VectorLayer that are covered
	 *         by the input MasonGeometry.
	 */
	public final List<MasonGeometry> touchingFeatures(MasonGeometry inputMasonGeometry) {
		List<MasonGeometry> touchingFeatures = new ArrayList<MasonGeometry>();
		Envelope envelope = inputMasonGeometry.getGeometry().getEnvelopeInternal();
		envelope.expandBy(Math.max(envelope.getHeight(), envelope.getWidth()) * 0.01);
		List<?> geometriesList = spatialIndex.query(envelope);

		if (inputMasonGeometry.preparedGeometry == null) {
			inputMasonGeometry.preparedGeometry = PreparedGeometryFactory.prepare(inputMasonGeometry.getGeometry());
		}

		for (final Object geometry : geometriesList) {
			final MasonGeometry otherMasonGeometry = (MasonGeometry) geometry;
			if (!inputMasonGeometry.equals(otherMasonGeometry)
					&& inputMasonGeometry.getGeometry().touches(otherMasonGeometry.getGeometry()))
				touchingFeatures.add(otherMasonGeometry);
		}
		return touchingFeatures;
	}

	/**
	 * Retrieves the geometries in the VectorLayer that contain the specified input
	 * Geometry.
	 *
	 * @param inputGeometry The Geometry to be contained.
	 * @return A List of MasonGeometry objects that contain the input Geometry in
	 *         the VectorLayer.
	 */
	public final List<MasonGeometry> containingFeatures(Geometry inputGeometry) {
		List<MasonGeometry> containingFeatures = new ArrayList<MasonGeometry>();
		Envelope envelope = inputGeometry.getEnvelopeInternal();
		envelope.expandBy(Math.max(envelope.getHeight(), envelope.getWidth()) * 0.01);
		List<?> geometriesList = spatialIndex.query(envelope);

		for (final Object geometry : geometriesList) {
			final MasonGeometry otherMasonGeometry = (MasonGeometry) geometry;
			Geometry otherGeometry = otherMasonGeometry.getGeometry();
			if (!inputGeometry.equals(otherGeometry) && otherGeometry.contains(inputGeometry))
				containingFeatures.add(otherMasonGeometry);
		}
		return containingFeatures;
	}

	/**
	 * Retrieves all the geometries in the VectorLayer that are fully contained
	 * within the boundaries of a given Geometry.
	 *
	 * @param inputGeometry The reference geometry defining the containment
	 *                      boundary.
	 * @return A List containing the geometries that are entirely contained within
	 *         the specified input geometry.
	 */
	public final List<MasonGeometry> containedFeatures(Geometry inputGeometry) {
		final List<MasonGeometry> containedFeatures = new ArrayList<>();
		final Envelope envelope = inputGeometry.getEnvelopeInternal();
		envelope.expandBy(Math.max(envelope.getHeight(), envelope.getWidth()) * 0.01);
		final List<?> geometriesList = spatialIndex.query(envelope);

		for (final Object geometry : geometriesList) {
			final MasonGeometry otherMasonGeometry = (MasonGeometry) geometry;
			if (inputGeometry.contains(otherMasonGeometry.getGeometry()))
				containedFeatures.add(otherMasonGeometry);
		}
		return containedFeatures;
	}

	// Other

	/**
	 * Reads a Shapefile from the specified input directory and populates the
	 * provided VectorLayer with geometries. The Shapefile consists of two files: a
	 * .shp file containing geometry information and a .dbf file containing
	 * attribute data. The function reads both files and adds the geometries to the
	 * VectorLayer.
	 *
	 * @param urlShp The URL referring to the path to the .shp input file.
	 * @param urlDbf The URL referring to the path to the .dbf input file.
	 * @throws Exception If there is an error during Shapefile reading or geometry
	 *                   population.
	 */
	public static void readShapefile(URL urlShp, URL urlDbf, VectorLayer vectorLayer) throws Exception {
		ShapeFileImporter.read(urlShp, urlDbf, vectorLayer, MasonGeometry.class);
	}

	public static void readGPKG(URL urlPkg, VectorLayer vectorLayer) throws Exception {
		GeoPackageImporter.read(urlPkg, vectorLayer);
	}

	public Envelope clipEnvelope;
	DrawInfo2D myInfo;
	public AffineTransform worldToScreen;
	public org.locationtech.jts.geom.util.AffineTransformation jtsTransform;

	/**
	 * Updates the transform using the provided DrawInfo2D object.
	 *
	 * @param info The DrawInfo2D object containing the necessary information for
	 *             the update
	 */
	public void updateTransform(DrawInfo2D info) {
		// need to update the transform
		if (!info.equals(myInfo)) {
			myInfo = info;
			// compute the transform between world and screen coordinates, and
			// also construct a geom.util.AffineTransform for use in hit-testing
			// later
			Envelope myMBR = getMBR();

			worldToScreen = GeometryUtilities.worldToScreenTransform(myMBR, info);
			jtsTransform = GeometryUtilities.getPortrayalTransform(worldToScreen, this, info.draw);

			Point2D p1 = GeometryUtilities.screenToWorldPointTransform(worldToScreen, info.clip.x, info.clip.y);
			Point2D p2 = GeometryUtilities.screenToWorldPointTransform(worldToScreen, info.clip.x + info.clip.width,
					info.clip.y + info.clip.height);

			clipEnvelope = new Envelope(p1.getX(), p2.getX(), p1.getY(), p2.getY());
		}
	}

}