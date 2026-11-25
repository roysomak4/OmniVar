import java.util.ArrayList;
import java.util.HashMap;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Map;
import java.util.Random;

/**
 * Writes java.util.HashMap iteration-order vectors for src/util/java_hash.
 *
 * Usage: java HashMapOrderVectors.java > tests/golden/java_hashmap_vectors.txt
 * Each line: "S" or "I", TAB, keys in insertion order (comma-separated), TAB, iteration order
 * as indices into the key list (comma-separated).
 */
public class HashMapOrderVectors {
    static void emit(String kind, List<String> keys, List<?> typed) {
        Map<Object, Integer> map = new HashMap<>();
        for (int i = 0; i < typed.size(); i++) map.put(typed.get(i), i);
        StringBuilder order = new StringBuilder();
        for (Map.Entry<Object, Integer> e : map.entrySet()) {
            if (order.length() > 0) order.append(',');
            order.append(e.getValue());
        }
        System.out.println(kind + "\t" + String.join(",", keys) + "\t" + order);
    }

    public static void main(String[] args) {
        Random rnd = new Random(7);
        String[] chroms = {"chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11",
                "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22",
                "chrX", "chrY", "chrM", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15",
                "16", "17", "18", "19", "20", "21", "22", "X", "Y", "MT", "chrEBV", "chr1_KI270706v1_random"};
        // Chromosome subsets in random order
        for (int t = 0; t < 200; t++) {
            List<String> pool = new ArrayList<>(List.of(chroms));
            java.util.Collections.shuffle(pool, rnd);
            List<String> keys = pool.subList(0, 1 + rnd.nextInt(pool.size()));
            emit("S", keys, keys);
        }
        // Random alphanumeric strings, up to a few thousand keys (resizes)
        for (int t = 0; t < 60; t++) {
            int n = 1 + rnd.nextInt(t < 50 ? 100 : 3000);
            LinkedHashSet<String> set = new LinkedHashSet<>();
            while (set.size() < n) {
                int len = 1 + rnd.nextInt(12);
                StringBuilder sb = new StringBuilder();
                for (int i = 0; i < len; i++) sb.append("ACGTN0123456789_:-chrXYZ".charAt(rnd.nextInt(24)));
                set.add(sb.toString());
            }
            List<String> keys = new ArrayList<>(set);
            emit("S", keys, keys);
        }
        // Integer positions: dense windows, sparse genome positions, negatives, and strided
        // values that collide in small tables (bins of >= 9 before resize).
        for (int t = 0; t < 120; t++) {
            LinkedHashSet<Integer> set = new LinkedHashSet<>();
            int n = 1 + rnd.nextInt(t < 100 ? 300 : 4000);
            int mode = t % 4;
            int base = rnd.nextInt(200_000_000);
            while (set.size() < n) {
                int v;
                if (mode == 0) v = base + rnd.nextInt(n * 2 + 10);
                else if (mode == 1) v = rnd.nextInt(250_000_000);
                else if (mode == 2) v = rnd.nextInt(4 * n + 2000) - (2 * n + 1000);
                else v = base + 16 * rnd.nextInt(n * 3 + 10);
                set.add(v);
            }
            List<Integer> typed = new ArrayList<>(set);
            List<String> keys = new ArrayList<>();
            for (Integer v : typed) keys.add(v.toString());
            emit("I", keys, typed);
        }
    }
}
