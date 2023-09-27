package sim.util.geo;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Objects;
import java.util.Random;
import java.util.stream.Collectors;

import sim.util.Bag;

/**
 * A class of various utilities.
 *
 */
public class Utilities {

	/** It sorts a Map on the basis of its values, and on the method provided (descending)
	 *
	 * @param map the map;
	 * @param descending if true, values are sorted in a descending order, otherwise ascending;
	 */
	public static <K, V extends Comparable<? super V>> Map<K, V> sortByValue(Map<K, V> map, boolean descending) {
		if (descending) return map.entrySet()
				.stream()
				.sorted(Map.Entry.comparingByValue(Collections.reverseOrder()))
				.collect(Collectors.toMap(
						Map.Entry::getKey,
						Map.Entry::getValue,
						(e1, e2) -> e1,
						LinkedHashMap::new
						));

		else return map.entrySet()
				.stream()
				.sorted(Map.Entry.comparingByValue())
				.collect(Collectors.toMap(
						Map.Entry::getKey,
						Map.Entry::getValue,
						(e1, e2) -> e1,
						LinkedHashMap::new
						));
	}
	
	/**
	 * It filters the content of the map on the basis of the key values. Only entries whose keys are contained in the given Bag are kept.
	 * A new map is returned.
	 *
	 * @param map, the input map;
	 * @param filter, a Bag containing the desired keys;
	 */
	public static HashMap<MasonGeometry, Double> filterMap(HashMap<MasonGeometry, Double> map, Bag filter) {
		HashMap<MasonGeometry, Double> mapFiltered = new HashMap<MasonGeometry, Double> (map);
		ArrayList<MasonGeometry> result = new ArrayList<MasonGeometry>();
		for(MasonGeometry key : mapFiltered.keySet()) {if(filter.contains(key)) result.add(key);}
		mapFiltered.keySet().retainAll(result);
		return mapFiltered;
	}

	/**
	 * Given a certain value, it returns the key of the first entry whose values is equal to the input.
	 *
	 * @param map, the input map;
	 * @param value, the input value whose key is of interest;
	 */
	public static <K, V> K getKeyFromValue(Map<K, V> map, V value) {
		for (Entry<K, V> entry : map.entrySet()) {
			if (Objects.equals(value, entry.getValue())) return entry.getKey();
		}
		return null;
	}

	/**
	 * It returns a random value from a distribution with a given mean and a standard deviation.
	 * If a direction is provided, only values higher or lower than the mean are returned.
	 *
	 * @param mean the distribution's mean;
	 * @param sd the distribution's standard deviation;
	 * @param direction either "left", only values lower than the mean, or "right", higher than the mean, otherwise pass null;
	 */
	public static double fromDistribution(double mean, double sd, String direction) 	{
		Random random = new Random();
		double result = random.nextGaussian()*sd+mean;
		if (direction != null) {
			if ((direction.equals("left")) && (result > mean)) result = mean;
			if ((direction.equals("right")) && (result < mean)) result = mean;
		}
		if (result <= 0.00) result = mean;
		return result;
	}
	
    // Filter the map to keep entries where the value is >= thresholdValue
    public static <K> Map<K, Integer> filterMapByMinValue(Map<K, Integer> map, int minValue) {
        return map.entrySet().stream()
                .filter(entry -> entry.getValue() >= minValue)
                .collect(Collectors.toMap(Map.Entry::getKey, Map.Entry::getValue));
    }
	
    // Filter the map to keep entries below the threshold value
    public static <K, V extends Comparable<Double>> Map<K, Double> filterMapByPercentile(Map<K, Double> map, double percentile) {
        
    	double percentileValue = calculatePercentileThreshold(map, percentile);
    	return filterMapByThreshold(map, percentileValue);
    }

    // Filter the map to keep entries below the threshold value
    private static <K, V extends Comparable<Double>> Map<K, Double> filterMapByThreshold(Map<K, Double> map, Double thresholdValue) {
        return map.entrySet().stream()
                .filter(entry -> entry.getValue().compareTo(thresholdValue) < 0)
                .collect(Collectors.toMap(Map.Entry::getKey, Map.Entry::getValue));
    }
    
    // Calculate the threshold value for a given percentile
    private static <K> double calculatePercentileThreshold(Map<K, Double> map, double percentile) {
        // Convert the map values to an array for percentile calculation
        double[] valuesArray = map.values().stream().mapToDouble(Double::doubleValue).toArray();

        // Sort the values array (ascending order)
        java.util.Arrays.sort(valuesArray);

        // Calculate the index corresponding to the specified percentile
        int index = (int) Math.ceil(percentile * valuesArray.length) - 1;

        // Return the threshold value
        return valuesArray[index];
    }
    
    public static <K, V> HashMap<K, V> filterMapByIndex(HashMap<K, V> map, ArrayList<K> keyList) {
    	Map<K, V> filteredMap = map.entrySet().stream()
		    .filter(entry -> keyList.contains(entry.getKey()))
		    .collect(Collectors.toMap(Map.Entry::getKey, Map.Entry::getValue));
    	return new HashMap<K, V> (filteredMap);
    }
    	
}





