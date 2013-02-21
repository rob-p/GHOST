package edu.umd.cbcb.mathtools;

import java.util.Arrays;

public class COOMatVecMult extends MatVecMult {

    private int[] rptr;
    private int[] cptr;
    private double[] dat;
    private int nrow;
    private int ncol;

    public COOMatVecMult( int[] _rptr, int[] _cptr, double[] _dat, int _nrow, int _ncol ) {
	rptr = _rptr;
	cptr = _cptr;
	dat = _dat;
	ncol = _ncol;
	nrow = _nrow;
	assert( rptr.length == cptr.length );
	assert( cptr.length == dat.length );
    }

    @Override public void mult( double[] x, int offset1, int len1, int offset2, int len2 ) {
	// Zero out the result portion of the array
	Arrays.fill(x, offset2, offset2+len2, 0.0);
	// Perfrom the spmv 
	int N = dat.length;
	for (int i = 0; i < N; i++) {
	    x[offset2 + rptr[i]] += dat[i] * x[ offset1 + cptr[i] ];
	}	    
    } 

    @Override public int numRows() { return nrow; }
}

