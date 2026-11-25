import com.astrazeneca.vardict.data.fishertest.FisherExact;

import java.util.Random;

/**
 * Writes test vectors for src/stats/fisher: VarDictJava's FisherExact (R fisher.test port) on
 * 2x2 tables. Columns: refFwd refRev altFwd altRev, then getPValue(), getOddRatio(),
 * getPValueLess(), getPValueGreater() (doubles as Double.toString).
 *
 * Usage: java -cp "tests/oracle/lib/*" tests/oracle/FisherVectors.java > tests/golden/fisher_vectors.tsv
 */
public class FisherVectors {
    static void emit(int a, int b, int c, int d) {
        String line;
        try {
            FisherExact f = new FisherExact(a, b, c, d);
            line = f.getPValue() + "\t" + f.getOddRatio() + "\t" + f.getPValueLess() + "\t" + f.getPValueGreater();
        } catch (Exception e) {
            line = "ERROR\t" + e.getClass().getSimpleName() + "\t\t";
        }
        System.out.println(a + "\t" + b + "\t" + c + "\t" + d + "\t" + line);
    }

    public static void main(String[] args) {
        for (int a = 0; a <= 8; a++)
            for (int b = 0; b <= 8; b++)
                for (int c = 0; c <= 8; c++)
                    for (int d = 0; d <= 8; d++)
                        if ((a + b + c + d) % 3 == 0 || a + b + c + d < 6) emit(a, b, c, d);
        Random rnd = new Random(11);
        for (int i = 0; i < 1500; i++) {
            int scale = new int[]{20, 100, 500, 3000}[i % 4];
            emit(rnd.nextInt(scale), rnd.nextInt(scale), rnd.nextInt(scale / 4 + 1), rnd.nextInt(scale / 4 + 1));
        }
        // typical variant tables: deep reference, few alt reads, strand bias
        for (int i = 0; i < 500; i++) {
            int ref = 50 + rnd.nextInt(2000);
            int alt = 1 + rnd.nextInt(60);
            int rf = rnd.nextInt(ref + 1);
            int af = rnd.nextInt(alt + 1);
            emit(rf, ref - rf, af, alt - af);
        }
    }
}
