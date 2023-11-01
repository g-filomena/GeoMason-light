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
package sim.graph;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.Map;

import org.locationtech.jts.geom.Geometry;

import pedSim.cognitiveMap.Region;
import pedSim.engine.PedSimCity;
import sim.field.geo.VectorLayer;
import sim.util.geo.AttributeValue;
import sim.util.geo.MasonGeometry;

/**
 * A class representing building metadata within a geographic information
 * system. Each instance of this class stores information about a building,
 * including its unique identifier, land use classification, geographic
 * geometry, associated node in a graph, local and global landmarkness scores,
 * and a label based on the Urban DMA categorisation.
 *
 * This class is used to organize and manage data related to buildings in a
 * geographical context, facilitating various spatial analysis and modeling
 * tasks.
 */
public class Building {

	public int buildingID; // Unique identifier for the building
	public String landUse; // Classification of land use for the building
	public MasonGeometry geometry; // Geographic geometry representing the building
	public NodeGraph node; // Associated node in a graph (if applicable)
	public String DMA; // Urban DMA categorisation

	public Map<String, AttributeValue> attributes = new HashMap<>();

	/**
	 * Returns all the buildings enclosed between two nodes.
	 *
	 * @param originNode      The first node.
	 * @param destinationNode The second node.
	 * @return A list of buildings.
	 */
	public static ArrayList<MasonGeometry> getBuildings(NodeGraph originNode, NodeGraph destinationNode) {

		Geometry smallestCircle = GraphUtils.enclosingCircleBetweenNodes(originNode, destinationNode);
		return PedSimCity.buildings.containedFeatures(smallestCircle);
	}

	/**
	 * Get buildings within a specified region.
	 *
	 * @param region The region for which buildings are to be retrieved.
	 * @return An ArrayList of MasonGeometry objects representing buildings within
	 *         the region.
	 */
	public static ArrayList<MasonGeometry> getBuildingsWithinRegion(Region region) {
		VectorLayer regionNetwork = region.regionNetwork;
		Geometry convexHull = regionNetwork.getConvexHull();
		return PedSimCity.buildings.containedFeatures(convexHull);
	}
}
