package sim.graph;
import sim.util.geo.MasonGeometry;

/**
 * A class to store buildings' metainformation
 *
 */

public class Building {

	public int buildingID;
	public String landUse;
	public MasonGeometry geometry;
	public NodeGraph node;
	public double localLandmarkness;
	public double globalLandmarkness;
	public String DMA;
}
