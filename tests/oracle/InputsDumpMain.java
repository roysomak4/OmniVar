import com.astrazeneca.vardict.CmdParser;
import com.astrazeneca.vardict.Configuration;
import com.astrazeneca.vardict.VarDictLauncher;
import com.astrazeneca.vardict.data.Reference;
import com.astrazeneca.vardict.data.ReferenceResource;
import com.astrazeneca.vardict.data.Region;
import com.astrazeneca.vardict.data.scopedata.GlobalReadOnlyScope;
import com.astrazeneca.vardict.data.scopedata.InitialData;
import com.astrazeneca.vardict.modules.RecordPreprocessor;
import htsjdk.samtools.SAMRecord;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.lang.reflect.Field;
import java.lang.reflect.Method;
import java.nio.charset.StandardCharsets;
import java.util.*;

/**
 * Dumps what VarDictJava 1.8.3 derives from its inputs before CIGAR parsing, for comparison with
 * the C++ port's --dev-dump-inputs (M2 parity):
 *
 *   SAMPLE   sample  samplem
 *   AMPLICON parameters|null
 *   REGION   i  j  chr  start  end  gene  insertStart  insertEnd       (segments, in order)
 *   REF      chr:start-end  positions  digest  seedKeys  digest       (ReferenceResource.getReference)
 *   READS    bamIndex  chr:start-end  returned  totalReads  duplicateReads  digest  (RecordPreprocessor)
 *
 * Digests are FNV-1a 64 (lower-case hex) over:
 *   reference: "pos:base;" for positions in ascending order
 *   seeds:     "key=p1,p2,...;" for keys in ascending order, positions in list order
 *   reads:     "name\tflags\tstart\tcigar\n" for each record returned, in order
 *
 * Usage: java -cp "lib/*" InputsDumpMain.java <VarDict args...>
 * A -G value ending in .csv is read as a thin-FASTA CSV (as OracleMain does).
 */
public class InputsDumpMain {
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

    static final class Fnv {
        long h = 0xcbf29ce484222325L;

        void add(String s) {
            for (byte b : s.getBytes(StandardCharsets.UTF_8)) {
                h ^= (b & 0xff);
                h *= 0x100000001b3L;
            }
        }

        String hex() {
            return Long.toHexString(h);
        }
    }

    static String str(Object o) {
        return o == null ? "null" : o.toString();
    }

    @SuppressWarnings("unchecked")
    public static void main(String[] args) throws Exception {
        Configuration conf = new CmdParser().parseParams(args);
        ReferenceResource ref = (conf.fasta != null && conf.fasta.endsWith(".csv"))
                ? new CsvReference(conf.fasta) : new ReferenceResource();
        VarDictLauncher launcher = new VarDictLauncher(ref);
        Method init = VarDictLauncher.class.getDeclaredMethod("initResources", Configuration.class);
        init.setAccessible(true);
        init.invoke(launcher, conf);
        Field segField = VarDictLauncher.class.getDeclaredField("segments");
        segField.setAccessible(true);
        List<List<Region>> segments = (List<List<Region>>) segField.get(launcher);

        GlobalReadOnlyScope scope = GlobalReadOnlyScope.instance();
        StringBuilder out = new StringBuilder();
        out.append("SAMPLE\t").append(str(scope.sample)).append('\t').append(str(scope.samplem)).append('\n');
        out.append("AMPLICON\t").append(str(scope.ampliconBasedCalling)).append('\n');
        int i = 0;
        for (List<Region> list : segments) {
            int j = 0;
            for (Region r : list) {
                out.append("REGION\t").append(i).append('\t').append(j).append('\t').append(r.chr).append('\t')
                        .append(r.start).append('\t').append(r.end).append('\t').append(r.gene).append('\t')
                        .append(r.insertStart).append('\t').append(r.insertEnd).append('\n');
                j++;
            }
            i++;
        }
        List<String> bamGroups = new ArrayList<>();
        bamGroups.add(conf.bam.getBam1());
        if (conf.bam.hasBam2()) bamGroups.add(conf.bam.getBam2());

        for (List<Region> list : segments) {
            for (Region region : list) {
                String where = region.chr + ":" + region.start + "-" + region.end;
                try {
                    Reference reference = ref.getReference(region);
                    Fnv refFnv = new Fnv();
                    for (Integer pos : new TreeSet<>(reference.referenceSequences.keySet())) {
                        refFnv.add(pos + ":" + reference.referenceSequences.get(pos) + ";");
                    }
                    Fnv seedFnv = new Fnv();
                    for (String key : new TreeSet<>(reference.seed.keySet())) {
                        StringBuilder sb = new StringBuilder(key).append('=');
                        List<Integer> ps = reference.seed.get(key);
                        for (int k = 0; k < ps.size(); k++) {
                            if (k > 0) sb.append(',');
                            sb.append(ps.get(k));
                        }
                        seedFnv.add(sb.append(';').toString());
                    }
                    out.append("REF\t").append(where).append('\t').append(reference.referenceSequences.size())
                            .append('\t').append(refFnv.hex()).append('\t').append(reference.seed.size())
                            .append('\t').append(seedFnv.hex()).append('\n');
                } catch (Exception e) {
                    Throwable t = e;
                    out.append("REF\t").append(where).append("\tERROR\t").append(t.getClass().getSimpleName())
                            .append('\n');
                }
                for (int b = 0; b < bamGroups.size(); b++) {
                    try {
                        RecordPreprocessor pre = new RecordPreprocessor(bamGroups.get(b).split(":"), region,
                                new InitialData());
                        Fnv fnv = new Fnv();
                        int returned = 0;
                        SAMRecord rec;
                        while ((rec = pre.nextRecord()) != null) {
                            fnv.add(rec.getReadName() + "\t" + rec.getFlags() + "\t" + rec.getAlignmentStart() + "\t"
                                    + rec.getCigarString() + "\n");
                            returned++;
                        }
                        pre.close();
                        out.append("READS\t").append(b + 1).append('\t').append(where).append('\t').append(returned)
                                .append('\t').append(pre.totalReads).append('\t').append(pre.duplicateReads)
                                .append('\t').append(fnv.hex()).append('\n');
                    } catch (Exception e) {
                        out.append("READS\t").append(b + 1).append('\t').append(where).append("\tERROR\t")
                                .append(e.getClass().getSimpleName()).append('\n');
                    }
                }
            }
        }
        System.out.print(out);
        System.out.flush();
        System.exit(0);
    }
}
