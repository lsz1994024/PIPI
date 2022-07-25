package proteomics.FM;

import java.util.Arrays;
        import java.util.Comparator;
        import java.util.List;

// code with some modification is from https://algosome.com/articles/sort-array-index-java.html

/**
 * Class to sort the indexes of an array based upon their values. Note the array or Collection passed
 * into the constructor is not itself sorted.
 * doubles,
 * @author G, Cope
 *
 */
public class IndexSorter implements Comparator<Integer>{

    private final char[] values;

    private final Integer[] indexes;

    /**
     * Constructs a new IndexSorter based upon the parameter array.
     * @param text
     */
    public IndexSorter(char[] text, Integer[] indexes){
        this.values = text;
        this.indexes = indexes;
    }


    /**
     * Sorts the underlying index array based upon the values provided in the constructor. The underlying value array is not sorted.
     */
    public void sort(){
        // This uses TimSort Algorithm which has O(n * log(n)) time complexity
        // https://docs.oracle.com/javase/7/docs/api/java/util/Collections.html
        // https://en.wikipedia.org/wiki/Timsort
        Arrays.sort(indexes, this);
    }

    /**
     * Compares the two values at index arg0 and arg0
     * @param arg0 The first index
     * @param arg1 The second index
     * @return The result of calling compareTo on T objects at position arg0 and arg1
     */
    @Override
    public int compare(Integer arg0, Integer arg1) {
        return Character.compare(values[arg0], values[arg1]);
    }

}

