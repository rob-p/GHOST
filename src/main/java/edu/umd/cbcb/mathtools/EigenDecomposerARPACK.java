package edu.umd.cbcb.mathtools;

import org.netlib.arpack.ARPACK;
import org.netlib.util.StringW;
import org.netlib.util.booleanW;
import org.netlib.util.doubleW;
import org.netlib.util.floatW;
import org.netlib.util.intW;

import org.apache.commons.math.linear.*;
import org.apache.commons.math.random.RandomDataImpl;

import java.lang.Math;
import java.util.Arrays;
import java.util.ArrayList;

public class EigenDecomposerARPACK {

    /**
     * Returns the requested eigendecomposition obtained via the application of
     * the supplied MatVecMult if successful.
     *
     * @param mat Provides a matrix vector product (including the
     * associated matrix) for the Arnoldi iterations of the
     * eigenvector computation.
     *
     * @param _k The number of eigenpairs to compute
     *
     * @param _which The portion of the spectrum for which eigenpairs should be computed
     *
     *
     * @return An EigenDecompositionResult object containing the requested eigenpairs.
     *
     */
    public static EigenDecompositionResult decompose( MatVecMult mat, int _k, String _which ) {
	ARPACK apack = ARPACK.getInstance();
	String[] errors = {"Normal exit.", 
			   "N must be positive.",
			   "NEV must be positive.", 
			   "NCV must be greater than NEV and less than or equal to N.",
			   "WHICH must be one of 'LM', 'SM', 'LA', 'SA' or 'BE'.", 
			   "BMAT must be one of 'I' or 'G'.", 
			   "Length of private work WORKL array is not sufficient.",
			   "Error return from trid. eigenvalue calculation", 
			   "Information error from LAPACK routine dsteqr.", 
			   "Starting vector is zero.",
			   "IPARAM(7) must be 1,2,3,4,5.", 
			   "IPARAM(7) = 1 and BMAT = 'G' are incompatible.", 
			   "NEV and WHICH = 'BE' are incompatible.",
			   "DSAUPD did not find any eigenvalues to sufficient accuracy.", 
			   "HOWMNY must be one of 'A' or 'S' if RVEC = .true.", 
			   "HOWMNY = 'S' not yet implemented"};
 
	intW nev = new intW(_k);
	int nevv = nev.val;
	int lenAF = 2 * _k;
	int ncv = Math.min( mat.numRows(), lenAF );
	int maxncv = Math.min( mat.numRows(), ncv+2 );
	String bmat = "I";
	String which = _which;
    
	int n = mat.numRows();
	int maxn = n;
    
	double[] workl = new double[3*maxncv*maxncv + 6*maxncv];
	double[] d = new double[3*maxncv];
	double[] resid = new double[maxn];
	double[] v = new double[maxn*maxncv];
	double[] workd = new double[3*maxn];
      
	doubleW tol = new doubleW(1e-8);
	double told = tol.val;
	intW ido = new intW(0);
	intW info = new intW(0);
	intW ierr = new intW(0);
	int ishifts = 1;
	int maxitr = 300;
	int model = 1;
	int[] iparam = new int[11];
	int[] ipntr = new int[14];
	boolean[] select = new boolean[maxncv];
	int ldv = maxn;
	int lworkl = ncv * (ncv + 8);
    
	int arg21 = 1;
	intW arg22 = new intW(2);

	double sigma = 0.0;
	iparam[0] = ishifts;
	iparam[2] = maxitr;
	iparam[6] = model;
      
	int itn = 0;
	boolean iter = true;
	do {
	    apack.dsaupd(ido, bmat, n, which, nevv, tol, resid, ncv, v, ldv,
			 iparam, ipntr, workd, workl, lworkl, info); //, arg21, arg22);

	    if ( ido.val == -1 || ido.val == 1 ) {

		int s1 = ipntr[0]-1;
		int s2 = ipntr[1]-1;
		mat.mult( workd, s1, n, s2, n );

	    } else {
		iter = false;
	    }
	    itn += 1;
	} while ( iter );
      
	// check for convergence or error
	if (info.val < 0) {
	    System.out.println("ARPACK ERROR : "+errors[-info.val]);
	} else {
	    boolean rvec = true;
	    nev = new intW(nevv);
	    told = tol.val;
	    apack.dseupd(rvec, "All", select, d, v, ldv, sigma, bmat, n, which,
			 nev, told, resid, ncv, v, ldv, iparam, ipntr, workd,
			 workl, lworkl, ierr);//, arg21, arg22);
	}
      
	double[] ew = null;
	ArrayList<double[]> ev = null;

	boolean success = false;
	EigenDecompositionResult res = new EigenDecompositionResult(success, ew, ev);
	if (ierr.val != 0) {
	    System.out.println("Got back IERR = "+ierr.val);
	    success = false;
	} else {
	    int numEigenVec = nev.val;
	    int nconv = iparam[4];
	    ew = Arrays.copyOfRange(d, 0, nconv);
	    ev = new ArrayList<double[]>(nconv);
	    for (int i=0; i < nconv; ++i) {
		ev.add( Arrays.copyOfRange(v, n*i, n*(i+1) ) );
	    }
	    success = true;
	    res = new EigenDecompositionResult(success, ew, ev);
	}

	return res;
    }

    public static void main( String[] args ) { 
	RandomDataImpl rand = new RandomDataImpl();
	double density = Double.parseDouble(args[1]);
	int N = Integer.parseInt(args[0]);
	
	ArrayList<Integer> rowinds = new ArrayList<Integer>();
	ArrayList<Integer> colinds = new ArrayList<Integer>();
	ArrayList<Double> data = new ArrayList<Double>();
	
	OpenMapRealMatrix M = new OpenMapRealMatrix(N,N);
	for (int i = 0; i < N; ++i) {
	    for (int j = i; j < N; ++j) { 
		double nr = rand.nextUniform(0.0, 1.0);
		if (nr > (1.0 - density)) { 
		    double e = rand.nextUniform(0.1,1.0);
		    M.setEntry(i,j,e);
		    M.setEntry(j,i,e);
		    
		    rowinds.add(i); 
		    colinds.add(j);
		    data.add(e);
		    
		    if ( i != j ) {
			rowinds.add(j);
			colinds.add(i);
			data.add(e); 
		    }

		}
	    }
	}

	System.out.println("              Input Matrix                ");
	System.out.println("=========================================");
	System.out.println(M);
	System.out.print("\n\n");


	int nnz = data.size();
	int[] rinds = new int[rowinds.size()];
	for (int i = 0; i < nnz; ++i ) { rinds[i] = rowinds.get(i).intValue(); }
	int[] cinds = new int[colinds.size()];
	for (int i = 0; i < nnz; ++i ) { cinds[i] = colinds.get(i).intValue(); }
	double[] dat = new double[data.size()];
	for (int i = 0; i < nnz; ++i ) { dat[i] = data.get(i).doubleValue(); }

	MatVecMult mvm = new COOMatVecMult(rinds, cinds, dat, N, N);
	EigenDecompositionResult res = EigenDecomposerARPACK.decompose(mvm, 4, "LA");
	
	String sstring;
	if ( res.getSuccess() ) { sstring = "Succeeded"; } else { sstring = "Failed"; }

	System.out.println("      Results computed with ARPACK       ");
	System.out.println("=========================================");

	System.out.println("Decomposition : "+sstring);
	System.out.println("Eigenvalues : "+Arrays.toString(res.getEvals()));
	System.out.println("Eigenvectors : ");
	ArrayList<double[]> evecs = res.getEvecs();
	for (int i=0; i < evecs.size(); ++i ) {
	    System.out.println( Arrays.toString(evecs.get(i)) );
	}

	System.out.println("\n\n");

	mvm = new ApacheMatVecMult(M);
	res = EigenDecomposerARPACK.decompose(mvm, 4, "LA");

	if ( res.getSuccess() ) { sstring = "Succeeded"; } else { sstring = "Failed"; }

	System.out.println("      Results computed with ARPACK       ");
	System.out.println("=========================================");

	System.out.println("Decomposition : "+sstring);
	System.out.println("Eigenvalues : "+Arrays.toString(res.getEvals()));
	System.out.println("Eigenvectors : ");
	evecs = res.getEvecs();
	for (int i=0; i < evecs.size(); ++i ) {
	    System.out.println( Arrays.toString(evecs.get(i)) );
	}

	System.out.print("\n\n");

	EigenDecompositionImpl ed = new EigenDecompositionImpl(M, 1e-8);
	System.out.println("Results computed with Apache Commons Math");
	System.out.println("=========================================");

	System.out.println("Eigenvalues : "+Arrays.toString(ed.getRealEigenvalues()));
	System.out.println("Eigenvectors : ");
	System.out.println(ed.getV());
    }
  
}
