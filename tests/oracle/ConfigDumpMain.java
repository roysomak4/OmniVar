import com.astrazeneca.vardict.CmdParser;
import com.astrazeneca.vardict.Configuration;

import java.io.BufferedReader;
import java.io.InputStreamReader;
import java.lang.reflect.Field;
import java.lang.reflect.Modifier;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;

/**
 * Prints the Configuration that VarDictJava 1.8.3's CmdParser builds from command lines, one
 * "key=value" line per field, so the C++ parser can be compared field by field.
 *
 * Usage: java -cp "lib/*" ConfigDumpMain < cases.txt
 * Each stdin line is one command line, arguments separated by TAB (no shell quoting). For each
 * line the output is "### <line>", then the fields sorted by name, or "ERROR <class>: <message>".
 * Lines that would print the help page (no options, -H) are not supported: help() calls exit.
 */
public class ConfigDumpMain {
    static String show(Object v) {
        if (v == null) return "null";
        if (v instanceof String) return "\"" + ((String) v).replace("\\", "\\\\").replace("\t", "\\t") + "\"";
        if (v instanceof List) {
            List<String> parts = new ArrayList<>();
            for (Object o : (List<?>) v) parts.add(show(o));
            return "[" + String.join(",", parts) + "]";
        }
        return String.valueOf(v);
    }

    static void dump(String[] args) throws Exception {
        Configuration c = new CmdParser().parseParams(args);
        Map<String, String> out = new TreeMap<>();
        for (Field f : Configuration.class.getDeclaredFields()) {
            if (Modifier.isStatic(f.getModifiers())) continue;
            f.setAccessible(true);
            Object v = f.get(c);
            String name = f.getName();
            if (name.equals("exceptionCounter")) continue;
            if (name.equals("bam")) {
                Configuration.BamNames b = (Configuration.BamNames) v;
                out.put("bam.raw", show(b == null ? null : b.getBamRaw()));
                out.put("bam.bam1", show(b == null ? null : b.getBam1()));
                out.put("bam.bam2", show(b == null ? null : b.getBam2()));
                out.put("bam.bamX", show(b == null ? null : b.getBamX()));
                continue;
            }
            if (name.equals("bedRowFormat")) {
                for (String col : new String[]{"chrColumn", "startColumn", "endColumn", "thickStartColumn",
                        "thickEndColumn", "geneColumn"}) {
                    Field cf = v.getClass().getDeclaredField(col);
                    cf.setAccessible(true);
                    out.put("bedRowFormat." + col, show(cf.get(v)));
                }
                continue;
            }
            out.put(name, show(v));
        }
        for (Map.Entry<String, String> e : out.entrySet()) {
            System.out.println(e.getKey() + "=" + e.getValue());
        }
    }

    public static void main(String[] argv) throws Exception {
        BufferedReader in = new BufferedReader(new InputStreamReader(System.in));
        String line;
        while ((line = in.readLine()) != null) {
            if (line.isEmpty()) continue;
            System.out.println("### " + line);
            try {
                dump(line.split("\t", -1));
            } catch (Exception e) {
                Throwable t = e;
                System.out.println("ERROR " + t.getClass().getSimpleName() + ": " + t.getMessage());
            }
        }
    }
}
