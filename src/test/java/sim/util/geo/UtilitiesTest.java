package sim.util.geo;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertNull;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.Map;
import java.util.Random;
import org.junit.jupiter.api.DisplayName;
import org.junit.jupiter.api.Test;

class UtilitiesTest {

  private static Map<String, Double> scores() {
    Map<String, Double> map = new LinkedHashMap<>();
    map.put("b", 2.0);
    map.put("c", 3.0);
    map.put("a", 1.0);
    return map;
  }

  @Test
  @DisplayName("sortByValue() returns a map iterating in sorted order")
  void sortByValueOrdersAscendingAndDescending() {
    assertEquals(Arrays.asList("a", "b", "c"),
        new ArrayList<>(Utilities.sortByValue(scores(), false).keySet()));
    assertEquals(Arrays.asList("c", "b", "a"),
        new ArrayList<>(Utilities.sortByValue(scores(), true).keySet()));
  }

  @Test
  void getKeyFromValueFindsTheKeyOrNull() {
    assertEquals("c", Utilities.getKeyFromValue(scores(), 3.0));
    assertNull(Utilities.getKeyFromValue(scores(), 99.0));
  }

  @Test
  @DisplayName("fromDistribution() honours the requested side of the mean")
  void fromDistributionRespectsDirection() {
    Random random = new Random(1234L);
    for (int draw = 0; draw < 500; draw++) {
      assertTrue(Utilities.fromDistribution(100.0, 30.0, "left", random) <= 100.0);
    }
    for (int draw = 0; draw < 500; draw++) {
      assertTrue(Utilities.fromDistribution(100.0, 30.0, "right", random) >= 100.0);
    }
  }

  @Test
  @DisplayName("fromDistribution() never returns a non-positive draw")
  void fromDistributionStaysPositive() {
    Random random = new Random(99L);
    for (int draw = 0; draw < 500; draw++) {
      assertTrue(Utilities.fromDistribution(1.0, 1000.0, null, random) > 0.0);
    }
  }

  @Test
  @DisplayName("fromDistribution() is reproducible under a seeded generator")
  void fromDistributionIsReproducibleWithASeed() {
    Random first = new Random(7L);
    Random second = new Random(7L);
    for (int draw = 0; draw < 50; draw++) {
      assertEquals(Utilities.fromDistribution(50.0, 10.0, null, first),
          Utilities.fromDistribution(50.0, 10.0, null, second), 0.0);
    }
  }

  @Test
  void filterMapByMinValueKeepsValuesAtOrAboveTheThreshold() {
    Map<String, Integer> counts = new LinkedHashMap<>();
    counts.put("a", 1);
    counts.put("b", 5);
    counts.put("c", 10);

    Map<String, Integer> filtered = Utilities.filterMapByMinValue(counts, 5);
    assertEquals(2, filtered.size());
    assertTrue(filtered.containsKey("b"));
    assertTrue(filtered.containsKey("c"));
  }

  @Test
  @DisplayName("filterMapByPercentile() keeps the values below the percentile")
  void filterMapByPercentileKeepsValuesBelowThreshold() {
    Map<String, Double> map = new LinkedHashMap<>();
    for (int value = 1; value <= 10; value++) {
      map.put("k" + value, (double) value);
    }
    // The 50th percentile of 1..10 is 5.0; entries strictly below it are 1..4.
    assertEquals(4, Utilities.filterMapByPercentile(map, 0.5).size());
    assertEquals(9, Utilities.filterMapByPercentile(map, 1.0).size());
  }

  @Test
  void filterMapByIndexKeepsTheListedKeys() {
    HashMap<String, Integer> map = new HashMap<>();
    map.put("a", 1);
    map.put("b", 2);
    map.put("c", 3);

    HashMap<String, Integer> filtered =
        Utilities.filterMapByIndex(map, new ArrayList<>(Arrays.asList("a", "c")));
    assertEquals(2, filtered.size());
    assertEquals(Integer.valueOf(1), filtered.get("a"));
    assertEquals(Integer.valueOf(3), filtered.get("c"));
  }
}
