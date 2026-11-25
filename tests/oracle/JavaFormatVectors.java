import java.text.DecimalFormat;
import java.util.LinkedHashSet;
import java.util.Random;
import java.util.Set;

/**
 * Writes test vectors for src/util/java_format: for each double, the outputs of the Java
 * formatting calls VarDictJava uses (DecimalFormat "0" .. "0.0000000", Utils.getRoundedValueToPrint,
 * Utils.roundHalfEven, String.format("%.Nf"), Double.toString, Math.round).
 *
 * Usage: java JavaFormatVectors.java > tests/golden/java_format_vectors.tsv
 * Columns: raw bits (hex), then the outputs; see the header line.
 */
public class JavaFormatVectors {
    static String pattern(int n) {
        return n == 0 ? "0" : "0." + "0".repeat(n);
    }

    // Copy of Utils.getRoundedValueToPrint (vardictjavasrc/.../Utils.java:100)
    static String roundedToPrint(String pattern, double value) {
        return value == Math.round(value)
                ? new DecimalFormat("0").format(value)
                : new DecimalFormat(pattern).format(value).replaceAll("0+$", "");
    }

    public static void main(String[] args) {
        Set<Double> values = new LinkedHashSet<>();
        double[] special = {0.0, -0.0, 1.0, -1.0, 0.5, 1.5, 2.5, -2.5, 0.125, 0.375, 0.0625, 0.05, 0.15, 0.25, 0.35,
                0.45, 0.55, 0.65, 0.75, 0.85, 0.95, 0.005, 0.015, 0.025, 0.0005, 0.00005, 0.000049, -0.0001, -0.0004,
                -0.0006, 0.1, 0.2, 0.3, 1e-3, 9.999e-4, 1e7, 9999999.999, 1e-7, 1.0E23, 2e23, 5e-324, 4.9e-324,
                Double.MIN_NORMAL, Double.MAX_VALUE, 123456789.0, 1234567.0, 0.001, 0.0011, 100.0, 22.5, 37.4,
                0.99999, 0.999995, 9.5, 10.5, 99.95, 99.995, 1.005, 2.675, 1.0000005, 0.30000000000000004,
                Double.NaN, Double.POSITIVE_INFINITY, Double.NEGATIVE_INFINITY, 2.82879384806159E17, 1.0E-5,
                4.35, 4.45, 0.0035, 0.0045, 1e300, 1e-300, 3.0E-44};
        for (double d : special) values.add(d);
        // Typical VarDict ratios: counts / depth, means
        for (int n = 1; n <= 64; n++) {
            for (int k = 0; k <= n; k++) {
                values.add((double) k / n);
                values.add((double) k / n * 40.0);
            }
        }
        // x.yyy5 decimal ties at several scales
        for (int i = 0; i < 600; i++) {
            values.add(i / 600.0 + 0.00005);
            values.add(i / 60.0 + 0.0005);
            values.add(i / 6.0 + 0.005);
        }
        Random rnd = new Random(42);
        for (int i = 0; i < 1000; i++) {
            values.add(rnd.nextDouble());
            values.add(rnd.nextDouble() * 100);
            values.add(-rnd.nextDouble() * 10);
            double bits = Double.longBitsToDouble(rnd.nextLong());
            if (Math.abs(bits) < 1e12 || i % 50 == 0) values.add(bits);  // keep the file small
            values.add(Math.pow(10, rnd.nextInt(40) - 20) * rnd.nextDouble());
        }

        StringBuilder hdr = new StringBuilder("bits");
        for (int n = 0; n <= 7; n++) hdr.append("\tdf").append(n);
        for (int n = 1; n <= 5; n++) hdr.append("\tprint").append(n);
        for (int n = 0; n <= 7; n++) hdr.append("\trhe").append(n);
        for (int n = 0; n <= 4; n++) hdr.append("\tfmt").append(n);
        hdr.append("\ttostring\tmathround");
        System.out.println(hdr);

        for (double v : values) {
            StringBuilder sb = new StringBuilder(Long.toHexString(Double.doubleToRawLongBits(v)));
            for (int n = 0; n <= 7; n++) sb.append('\t').append(new DecimalFormat(pattern(n)).format(v));
            for (int n = 1; n <= 5; n++) sb.append('\t').append(roundedToPrint(pattern(n), v));
            for (int n = 0; n <= 7; n++) {
                String s;
                try {
                    s = Double.toString(Double.parseDouble(new DecimalFormat(pattern(n)).format(v)));
                } catch (NumberFormatException e) {
                    s = "NFE";
                }
                sb.append('\t').append(s);
            }
            for (int n = 0; n <= 4; n++) sb.append('\t').append(String.format("%." + n + "f", v));
            sb.append('\t').append(Double.toString(v));
            sb.append('\t').append(Math.round(v));
            System.out.println(sb);
        }
    }
}
