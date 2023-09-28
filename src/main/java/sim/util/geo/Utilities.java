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

    /**
     * Sorts a Map based on its values.
     *
     * @param map        The map to be sorted.
     * @param descending If true, values are sorted in descending order; otherwise, they are sorted in ascending order.
     * @param <K>        The type of keys in the map.
     * @param <V>        The type of values in the map.
     * @return A sorted map based on values.
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
	
//    /**
//     * Filters a map based on the keys provided in a Bag. Only entries with keys in the Bag are kept.
//     *
//     * @param map    The input map.
//     * @param filter A Bag containing the desired keys.
//     * @return A new map containing only the filtered entries.
//     */
//	public static HashMap<MasonGeometry, Double> filterMap(HashMap<MasonGeometry, Double> map, Bag filter) {
//		HashMap<MasonGeometry, Double> mapFiltered = new HashMap<MasonGeometry, Double> (map);
//		ArrayList<MasonGeometry> result = new ArrayList<MasonGeometry>();
//		for(MasonGeometry key : mapFiltered.keySet()) 
//			if(filter.contains(key)) result.add(key);
//		mapFiltered.keySet().retainAll(result);
//		return mapFiltered;
//	}

    /**
     * Retrieves the key from a map based on a specified value.
     *
     * @param map   The input map.
     * @param value The value whose corresponding key is of interest.
     * @param <K>   The type of keys in the map.
     * @param <V>   The type of values in the map.
     * @return The key corresponding to the provided value.
     */
	public static <K, V> K getKeyFromValue(Map<K, V> map, V value) {
		for (Entry<K, V> entry : map.entrySet()) {
			if (Objects.equals(value, entry.getValue())) return entry.getKey();
		}
		return null;
	}

    /**
     * Generates a random value from a distribution with a given mean and standard deviation.
     * Optionally, values higher or lower than the mean can be returned based on the provided direction.
     *
     * @param mean      The mean of the distribution.
     * @param sd        The standard deviation of the distribution.
     * @param direction Either "left" (for values lower than the mean), "right" (for values higher than the mean), or null.
     * @return A random value from the distribution.
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
	
    /**
     * Filters a map to keep entries where the value is greater than or equal to the provided threshold value.
     *
     * @param map        The input map.
     * @param minValue   The minimum value for filtering.
     * @param <K>        The type of keys in the map.
     * @return A filtered map.
     */
    public static <K> Map<K, Integer> filterMapByMinValue(Map<K, Integer> map, int minValue) {
        return map.entrySet().stream()
                .filter(entry -> entry.getValue() >= minValue)
                .collect(Collectors.toMap(Map.Entry::getKey, Map.Entry::getValue));
    }
	
    /**
     * Filters a map to keep entries below the provided percentile threshold value.
     *
     * @param map       The input map.
     * @param percentile The percentile value (e.g., 0.95 for 95th percentile).
     * @param <K>       The type of keys in the map.
     * @return A filtered map.
     */
    public static <K, V extends Comparable<Double>> Map<K, Double> filterMapByPercentile(Map<K, Double> map, 
    		double percentile) {
    	double percentileValue = calculatePercentileThreshold(map, percentile);
    	return filterMapByThreshold(map, percentileValue);
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
    
    // Filter the map to keep entries below the threshold value
    private static <K, V extends Comparable<Double>> Map<K, Double> filterMapByThreshold(Map<K, Double> map, 
    		double thresholdValue) {
        return map.entrySet().stream()
                .filter(entry -> entry.getValue().compareTo(thresholdValue) < 0)
                .collect(Collectors.toMap(Map.Entry::getKey, Map.Entry::getValue));
    }
    
    /**
     * Filters a map to keep entries that match the provided key list.
     *
     * @param map     The input map.
     * @param keyList The list of keys to include in the filtered map.
     * @param <K>     The type of keys in the map.
     * @param <V>     The type of values in the map.
     * @return A filtered map.
     */
    public static <K, V> HashMap<K, V> filterMapByIndex(HashMap<K, V> map, ArrayList<K> keyList) {
    	Map<K, V> filteredMap = map.entrySet().stream()
		    .filter(entry -> keyList.contains(entry.getKey()))
		    .collect(Collectors.toMap(Map.Entry::getKey, Map.Entry::getValue));
    	return new HashMap<K, V> (filteredMap);
    }

}





