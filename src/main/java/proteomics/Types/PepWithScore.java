package proteomics.Types;

import java.util.Set;

public class PepWithScore implements Comparable<PepWithScore> {
    public String pepSeq;
    public double score;
    public double count;
    public boolean isDecoy = false;
    public boolean hasPTM;
    public String proteins;
    public PepWithScore(String pepSeq, double score, boolean isDecoy, boolean hasPTM, String proteins) {
        this.pepSeq = pepSeq;
        this.score = score;
        this.count = count;
        this.isDecoy = isDecoy;
        this.hasPTM = hasPTM;
        this.proteins = proteins;
    }

    public int compareTo(PepWithScore other) {
        if (score < other.score) {  //decreasing order
            return 1;
        } else {
            return -1;
        }
    }

}
