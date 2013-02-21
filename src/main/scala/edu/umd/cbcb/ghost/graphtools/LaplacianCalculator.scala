package edu.umd.cbcb.ghost.graphtools

import scala.math._
import scala.collection.JavaConversions._

import org.jgrapht.graph._
//import org.apache.mahout.math._
import no.uib.cipr.matrix.sparse.{ FlexCompRowMatrix => SparseRowMat }
import no.uib.cipr.matrix.DenseVector

import edu.umd.cbcb.ghost.io.{ SimpleVertex }

class LaplacianCalculator[ T <: AbstractGraph[SimpleVertex, SequenceWeightedEdge] ](G: T, w:(Double, Double) => Double) {

  def compute : (SparseRowMat, DenseVector) = {

    val N = G.vertexSet.size
    //println("Constructing Laplacian of size "+N+" x "+N)
    val diagWeights = new DenseVector(N)
    val wmat = new SparseRowMat( N, N ) //Array(N,N), true )
    val laplacian = new SparseRowMat( N, N )//Array(N,N), true )

    val vertIDMap = G.vertexSet.zipWithIndex.map{ vi => ( (vi._1, vi._2) )}.toMap

    // First, we build the weight matrix W and the diagonal D
    // Iterate over each vertex (row)
    G.vertexSet.foreach{
      vi =>
	var dweight = 0.0
        val cidx = vertIDMap(vi)
	// Iterate over v_i's adjacent edges
	G.edgesOf(vi).foreach{
	  eij =>
	    val (src, tgt, weight) = (G.getEdgeSource(eij), G.getEdgeTarget(eij), G.getEdgeWeight(eij))
	    val (i,j) = (vertIDMap(src), vertIDMap(tgt))
	    val pweight = w(weight, eij.seqWeight)
	    dweight += pweight

	    if (i == cidx) {
	      wmat.set(i, j, pweight)
	    } else {
	      wmat.set(j, i, pweight)
	    }
	}

      diagWeights.set(cidx, dweight)
    }

    // Now, we can form L = I - D^{-1/2} * W * D^{-1/2}

    // Iterate over each vertex (row)
    G.vertexSet.foreach{
      vi =>
	val cidx = vertIDMap(vi)
	var rowSum = 0.0
	// Iterate over v_i's adjacent edges
	G.edgesOf(vi).foreach{
	  eij =>
	    val (src, tgt) = (G.getEdgeSource(eij), G.getEdgeTarget(eij))
	    val (i, j) = (vertIDMap(src), vertIDMap(tgt))
	    val (di, dj) = (diagWeights.get(i), diagWeights.get(j) )
	    var ow = 0.0
	    if (i == cidx) {
	      ow = wmat.get(i,j) / (sqrt(di * dj))
	      laplacian.set(i, j, -ow)
	    } else {
	      ow = wmat.get(j,i) / (sqrt(di * dj))
	      laplacian.set(j, i, -ow)
	    }
	    rowSum  += ow
	}
      laplacian.set(cidx, cidx, 1.0)
    }

    return ( laplacian, diagWeights )
  }

}
