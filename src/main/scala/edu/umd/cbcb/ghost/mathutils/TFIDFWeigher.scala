package edu.umd.cbcb.ghost.mathutils

import scala.collection.parallel.{ ParIterable, ParMap => PMap }
import scala.collection.parallel.mutable.{ ParHashSet => PMHashSet, ParHashMap => PMHashMap, ParArray }
import scala.math.log

import org.ejml.simple.SimpleMatrix
import org.ejml.ops.NormOps
import org.ejml.data.DenseMatrix64F

object TFIDFWeighter{

  def computeIDFVector( lm:SimpleMatrix, rm:SimpleMatrix ) = {
    val numDocs = (lm.numRows + rm.numRows).toDouble

    val idf = Array.ofDim[Double]( lm.numCols )
    var inds = Range(0, idf.size).toSet.par

    def l0Norm( v : DenseMatrix64F ) = {
      var (s, size) = (0, v.getNumElements)
      var n = 0.0
      while ( s < size ){
	n += { if ( v.get(s) > 0 ) { 1.0 } else { 0.0 } }; s += 1
      }
      n
    }
  
    // compute the idf vector
    inds.foreach{
      i =>
	idf(i) = log( numDocs /
		     (1.0 + l0Norm( lm.extractVector(false,i).getMatrix )
		      + l0Norm( rm.extractVector(false,i).getMatrix ) )
		     
		     // (1.0 + NormOps.normP1( lm.extractVector(false,i).getMatrix )
		     //  + NormOps.normP1( rm.extractVector(false,i).getMatrix ) )
		   )
    }
    
    new SimpleMatrix( Array( idf.seq.toArray ) )
  }
  

}
