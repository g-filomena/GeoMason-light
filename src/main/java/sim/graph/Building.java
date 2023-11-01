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
import java.util.HashMap;
import java.util.Map;

import sim.util.geo.AttributeValue;
import sim.util.geo.MasonGeometry;

/**
 * A class representing building metadata within a geographic information system.
 * Each instance of this class stores information about a building, including its
 * unique identifier, land use classification, geographic geometry, associated node
 * in a graph, local and global landmarkness scores, and a label based on the Urban DMA categorisation.
 *
 * This class is used to organize and manage data related to buildings in a
 * geographical context, facilitating various spatial analysis and modeling tasks.
 */
public class Building {
    
	public int buildingID;           // Unique identifier for the building
    public String landUse;           // Classification of land use for the building
    public MasonGeometry geometry;   // Geographic geometry representing the building
    public NodeGraph node;           // Associated node in a graph (if applicable)
    public String DMA;               //  Urban DMA categorisation
    
	public Map<String, AttributeValue> attributes = new HashMap<>();
}

