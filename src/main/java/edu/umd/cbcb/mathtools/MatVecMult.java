package edu.umd.cbcb.mathtools;

public abstract class MatVecMult {
    public abstract void mult( double[] x, int offset1, int len1, int offset2, int len2 );
    public abstract int numRows();
}

