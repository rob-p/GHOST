package edu.umd.cbcb.mathtools;

import org.apache.commons.math.linear.*;

public class ApacheMatVecMult<MatT extends RealMatrix> extends MatVecMult {
    private MatT mat;
    public ApacheMatVecMult(MatT _mat) {
	mat = _mat;
    }

    @Override public void mult( double[] x, int offset1, int len1, int offset2, int len2 ) {    
	ArrayRealVector v1 = new ArrayRealVector(x, offset1, len1);
	double[] d = mat.operate(v1.getDataRef());
	for (int i = 0; i < len2; ++i) { x[offset2+i] = d[i]; }
    }
	
    @Override public int numRows() { return mat.getRowDimension(); } 

    public MatT getMat() { return mat; }
}
