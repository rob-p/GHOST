package edu.umd.cbcb.mathtools;

import java.util.ArrayList;

public class EigenDecompositionResult {

    private boolean success;
    private double[] evals;
    private ArrayList<double[]> evecs;
    
    public EigenDecompositionResult( boolean _success, double[] _evals, ArrayList<double[]> _evecs ) {
	success = _success;
	evals = _evals;
	evecs = _evecs;
    }
	
    public double[] getEvals() { return evals; }
    public ArrayList<double[]> getEvecs() { return evecs; }
    public boolean getSuccess() { return success; }
}
