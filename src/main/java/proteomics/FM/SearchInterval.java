package proteomics.FM;

public class SearchInterval {
    public int sp;
    public int ep;
    public boolean settled = true;
    public int matchedPos = 0;
    public SearchInterval(int sp, int ep){
        this.sp = sp;
        this.ep = ep;
    }

}