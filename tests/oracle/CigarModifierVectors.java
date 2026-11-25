import com.astrazeneca.vardict.Configuration;
import com.astrazeneca.vardict.data.ModifiedCigar;
import com.astrazeneca.vardict.data.Reference;
import com.astrazeneca.vardict.data.Region;
import com.astrazeneca.vardict.data.scopedata.GlobalReadOnlyScope;

import java.io.PrintStream;
import java.lang.reflect.Constructor;
import java.lang.reflect.Method;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Random;

/**
 * Writes test vectors for src/core/cigar_modifier: the jar's CigarModifier.modifyCigar (local
 * realignment, -k 1) on random reads. Each read is built from a random reference along a random
 * CIGAR (soft clips that continue the reference or not, end indels, close indel clusters,
 * mismatches near the ends, clips that are the reverse complement of nearby reference), so every
 * rewrite rule of modifyCigar fires many times. The reference may stop short of the read.
 *
 * Columns: position cigar query quality refStart refSeq indel maxReadLength, then the result:
 * position cigar query quality, and the number of exceptions modifyCigar reported (0/1).
 *
 * Usage: java -cp "tests/oracle/lib/*" tests/oracle/CigarModifierVectors.java 10000 > tests/golden/cigar_modifier_vectors.tsv
 * (the checked-in file has 10,000 cases; 60,000 were also checked once).
 */
public class CigarModifierVectors {
    static final char[] NT = {'A', 'C', 'G', 'T'};
    static Random rnd = new Random(20261101);

    static char comp(char c) {
        switch (c) {
            case 'A': return 'T';
            case 'T': return 'A';
            case 'C': return 'G';
            case 'G': return 'C';
            default: return c;
        }
    }

    static int len(int small, int big) {
        return rnd.nextInt(3) == 0 ? 1 + rnd.nextInt(big) : 1 + rnd.nextInt(small);
    }

    public static void main(String[] args) throws Exception {
        Configuration config = new Configuration();
        GlobalReadOnlyScope.init(config, new HashMap<>(), "", "", null, new HashMap<>(), new HashMap<>());
        Constructor<?> ctor = Class.forName("com.astrazeneca.vardict.modules.CigarModifier").getDeclaredConstructor(
                int.class, String.class, String.class, String.class, Reference.class, int.class, Region.class, int.class);
        ctor.setAccessible(true);
        Method modify = ctor.getDeclaringClass().getMethod("modifyCigar");

        PrintStream out = System.out;
        PrintStream quiet = new PrintStream(new java.io.OutputStream() {
            public void write(int b) {}
        });
        int n = args.length > 0 ? Integer.parseInt(args[0]) : 10000;
        for (int c = 0; c < n; c++) {
            // reference
            int refStart = 1 + rnd.nextInt(60);
            int refLen = 340 - rnd.nextInt(60) - refStart;
            StringBuilder refSeq = new StringBuilder();
            for (int i = 0; i < refLen; i++) {
                // some low-complexity stretches
                if (i > 0 && rnd.nextInt(6) == 0) refSeq.append(refSeq.charAt(i - 1));
                else refSeq.append(NT[rnd.nextInt(4)]);
            }
            Reference ref = new Reference();
            for (int i = 0; i < refLen; i++) {
                ref.referenceSequences.put(refStart + i, refSeq.charAt(i));
            }
            for (int i = 0; i + 12 <= refLen; i++) {
                ref.seed.computeIfAbsent(refSeq.substring(i, i + 12), k -> new ArrayList<>()).add(refStart + i);
            }
            int position = 70 + rnd.nextInt(120);

            // CIGAR
            List<int[]> cig = new ArrayList<>();  // {len, op char}
            int lead = rnd.nextInt(10);
            if (lead < 3) cig.add(new int[]{len(8, 35), 'S'});
            else if (lead == 3) cig.add(new int[]{len(3, 6), 'D'});
            else if (lead == 4) cig.add(new int[]{len(3, 6), 'I'});
            if (lead < 3 && rnd.nextInt(3) == 0) cig.add(new int[]{len(3, 6), rnd.nextBoolean() ? 'I' : 'D'});
            int segs = 1 + rnd.nextInt(5);
            for (int s = 0; s < segs; s++) {
                cig.add(new int[]{rnd.nextInt(2) == 0 ? len(9, 16) : 10 + rnd.nextInt(60), 'M'});
                if (s + 1 < segs) {
                    int k = rnd.nextInt(10);
                    if (k < 4) cig.add(new int[]{len(3, 8), 'D'});
                    else if (k < 8) cig.add(new int[]{len(3, 8), 'I'});
                    else if (k == 8) cig.add(new int[]{len(3, 8), 'N'});
                    else {  // two indels in a row
                        cig.add(new int[]{len(3, 4), 'D'});
                        cig.add(new int[]{len(3, 4), rnd.nextBoolean() ? 'I' : 'D'});
                    }
                }
            }
            int tail = rnd.nextInt(10);
            if (tail < 2) cig.add(new int[]{len(3, 6), rnd.nextBoolean() ? 'I' : 'D'});
            if (tail >= 2 && tail < 5 || tail == 0 && rnd.nextBoolean()) cig.add(new int[]{len(8, 35), 'S'});

            // read along the CIGAR
            StringBuilder q = new StringBuilder();
            int rp = position;
            boolean revcompClip = rnd.nextInt(8) == 0;
            for (int e = 0; e < cig.size(); e++) {
                int l = cig.get(e)[0];
                char op = (char) cig.get(e)[1];
                if (op == 'M') {
                    double mm = rnd.nextInt(4) == 0 ? 0.15 : 0.02;
                    for (int i = 0; i < l; i++) {
                        Character r = ref.referenceSequences.get(rp);
                        char b = r == null || rnd.nextDouble() < mm ? NT[rnd.nextInt(4)] : r;
                        q.append(b);
                        rp++;
                    }
                } else if (op == 'I') {
                    for (int i = 0; i < l; i++) q.append(NT[rnd.nextInt(4)]);
                } else if (op == 'D' || op == 'N') {
                    rp += l;
                } else if (op == 'S') {
                    boolean first = e == 0;
                    int mode = rnd.nextInt(3);
                    for (int i = 0; i < l; i++) {
                        char b;
                        if (revcompClip && first) {
                            Character r = ref.referenceSequences.get(position + 40 + l - 1 - i);
                            b = r == null ? 'N' : comp(r);
                        } else if (mode == 0) {
                            Character r = ref.referenceSequences.get(first ? position - l + i : rp + i);
                            b = r == null ? 'A' : r;
                            if (rnd.nextInt(10) == 0) b = NT[rnd.nextInt(4)];
                        } else if (mode == 1 && i < 3) {
                            Character r = ref.referenceSequences.get(first ? position - l + i : rp + i);
                            b = r == null ? 'C' : r;
                        } else {
                            b = NT[rnd.nextInt(4)];
                        }
                        q.append(b);
                    }
                }
            }
            StringBuilder qual = new StringBuilder();
            for (int i = 0; i < q.length(); i++) {
                qual.append((char) (33 + (rnd.nextInt(8) == 0 ? rnd.nextInt(12) : 20 + rnd.nextInt(21))));
            }
            // occasionally a quality string that is too short (exceptions inside modifyCigar)
            if (rnd.nextInt(50) == 0 && qual.length() > 3) qual.setLength(qual.length() - 3);
            StringBuilder cigar = new StringBuilder();
            int indel = 0;
            for (int[] el : cig) {
                cigar.append(el[0]).append((char) el[1]);
                if (el[1] == 'I' || el[1] == 'D') indel += el[0];
            }
            int maxReadLength = rnd.nextBoolean() ? 150 : q.length();

            Region region = new Region("chr1", 1, 500, "");
            config.exceptionCounter.set(0);
            System.setErr(quiet);
            ModifiedCigar r = (ModifiedCigar) modify.invoke(ctor.newInstance(position, cigar.toString(), q.toString(),
                    qual.toString(), ref, indel, region, maxReadLength));
            System.setErr(System.err);
            out.println(position + "\t" + cigar + "\t" + q + "\t" + qual + "\t" + refStart + "\t" + refSeq + "\t" + indel
                    + "\t" + maxReadLength + "\t" + r.position + "\t" + r.cigar + "\t" + r.querySequence + "\t"
                    + r.queryQuality + "\t" + config.exceptionCounter.get());
        }
    }
}
