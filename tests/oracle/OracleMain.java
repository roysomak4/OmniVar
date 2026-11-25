import com.astrazeneca.vardict.CmdParser;
import com.astrazeneca.vardict.Configuration;
import com.astrazeneca.vardict.VarDictLauncher;
import com.astrazeneca.vardict.data.ReferenceResource;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.HashMap;
import java.util.Map;
import java.util.NavigableMap;
import java.util.TreeMap;

/**
 * Runs VarDictJava 1.8.3 with a reference read from a "thin FASTA" CSV (chr,start,end,seq),
 * reproducing the integration-test harness without Mockito.
 *
 * Usage: java -cp "lib/*:." OracleMain <ref.csv> <VarDict args...>
 * The value given to -G is ignored for lookups; pass the CSV path there too.
 */
public class OracleMain {
    static class CsvReference extends ReferenceResource {
        private final Map<String, NavigableMap<Integer, String>> slices = new HashMap<>();
        private final String path;

        CsvReference(String path) throws IOException {
            this.path = path;
            try (BufferedReader in = new BufferedReader(new FileReader(path))) {
                String line;
                while ((line = in.readLine()) != null) {
                    if (line.isEmpty()) continue;
                    String[] f = line.split(",");
                    slices.computeIfAbsent(f[0], k -> new TreeMap<>()).put(Integer.parseInt(f[1]), f[3]);
                }
            }
        }

        @Override
        public String[] retrieveSubSeq(String fasta, String chr, int start, int end) {
            String region = ">" + chr + ":" + start + "-" + end;
            NavigableMap<Integer, String> byStart = slices.get(chr);
            if (byStart != null) {
                Map.Entry<Integer, String> e = byStart.floorEntry(start);
                if (e != null) {
                    int csvStart = e.getKey();
                    String seq = e.getValue();
                    int csvEnd = csvStart + seq.length() - 1;
                    if (csvEnd >= end) {
                        return new String[]{region, seq.substring(start - csvStart, end - csvStart + 1)};
                    }
                }
            }
            throw new RuntimeException("No region: " + region + " in file " + path);
        }
    }

    public static void main(String[] args) throws Exception {
        String csv = args[0];
        String[] vardictArgs = new String[args.length - 1];
        System.arraycopy(args, 1, vardictArgs, 0, vardictArgs.length);
        Configuration conf = new CmdParser().parseParams(vardictArgs);
        new VarDictLauncher(new CsvReference(csv)).start(conf);
    }
}
