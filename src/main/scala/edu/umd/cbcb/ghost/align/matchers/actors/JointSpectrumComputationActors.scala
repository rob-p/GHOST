package edu.umd.cbcb.ghost.matchers.actors

import scala.collection.mutable.{ HashSet => MHashSet, HashMap => MHashMap, MultiMap }
import scala.collection.immutable.{ TreeMap }
import scala.collection.JavaConversions._
import scala.actors.Actor
import scala.actors._
import scala.math._

import org.jgrapht.graph._
import org.jgrapht.traverse._
import org.ejml.simple._
import org.ejml.alg.dense.decomposition._
import org.apache.mahout.math._

import java.io.{ File => JFile, FileWriter => JFileWriter }
import java.io.{ FileOutputStream => JFileOutputStream, FileInputStream => JFileInputStream }
import java.io.{ BufferedOutputStream => JBufferedOutputStream }
import java.io.{ DataOutputStream => JDataOutputStream, DataInputStream => JDataInputStream }
import java.util.concurrent.{ ConcurrentHashMap => PHMap }

import net.robpatro.utils.console.ProgressBar
import edu.umd.cbcb.ghost.align.{ SignatureMap }
import edu.umd.cbcb.ghost.io.{ SimpleVertex }
import edu.umd.cbcb.ghost.graphtools.{ LaplacianCalculator, SequenceWeightedEdge }
import edu.umd.cbcb.ghost.mathtools.{ Histogram }

/** Messages accepted by SubgraphExtractor and SubgraphExtractorReceiver **/

// Tells an extractor to extract the subgraph centered at v and return the descriptor
case class ComputeJointSpectrum( u: String, v : String )

case class JointSpectrumMsg( u: String, v: String, s: Array[Double] )

// Tells the actor that they are done
case class Stop()

class JointSpectrumReceiver( val reqNumSpec : Int ) extends Actor {

  var M = Array.ofDim[Double](reqNumSpec, 100)
  var numReceivedSpec = 0
  val pb = new ProgressBar( reqNumSpec, "*" )
  val specMap = new MHashMap[(String, String), SimpleMatrix]

  def act() = {

    loopWhile( !done ) {
      react {
	case JointSpectrumMsg( u: String, v: String, s: Array[Double] ) => {
	  val sv = new SimpleMatrix(Array(s))
	  if (u < v) {
	    specMap( (u,v) ) = sv; //new SimpleMatrix( Array(s) )
	  } else {
	    specMap( (v,u) ) = sv; //new SimpleMatrix( Array(s) )
	  }
	  pb.update(specMap.size)
	}
	case Stop() => {
	  pb.done()
	  exit()
	}
      } // end react
    } // end loopWhile

  }

  def done() = { specMap.size == reqNumSpec }

}


class JointSpectrumComputer( val S: SignatureMap,  val level: Int, val rec : PHMap[(String,String), SimpleMatrix])  extends Actor {
  /**
   * The graph and k hop neighborhood with which
   * this actor will work
   **/

  def act() = {
    loop {
      react {
	case ComputeJointSpectrum( u, v ) => {
	  val spec = new SimpleMatrix(Array(generateJointSpectrum(u,v)))
	  if (u < v) {
	    rec += ( (u,v) -> spec )
	  } else {
	    rec += ( (v,u) -> spec )
	  }
	  //rec ! JointSpectrumMsg( u, v, generateJointSpectrum( u, v ) )
	}
	case Stop => exit()
      }
    }
  }


  def generateJointSpectrum( u: String, v: String ): Array[Double] = {
    type SubgraphT = S.SubgraphT

    // Get the neighborhoods of the constitutent vertices at the desired levels
    val neighborhoodU = S.neighborhoodByName(u, level)
    val neighborhoodV = S.neighborhoodByName(v, level)

    // Compute the joint subgraph
    val subGraph = new SubgraphT(S.G, neighborhoodU & neighborhoodV)
    val N = subGraph.vertexSet.size
    if ( N > 0) {
      // Compute the laplacian matrix of this subgraph
      val lapCalc = new LaplacianCalculator(subGraph, (x, y) => x )
      val (lap, vols) = lapCalc.compute

      /** Use the dense eigensolver **/
      import org.ejml.ops._
      import org.ejml.data.DenseMatrix64F
      import org.ejml.ops.EigenOps

      val denseL = new DenseMatrix64F(N, N)
      for (i <- 0 until N;
	   j <- i until N) {
	denseL.unsafe_set(i, j, lap.get(i,j))
	denseL.unsafe_set(j, i, lap.get(j,i))
      }

      // compute the eigendecomposition
      val decomposer = DecompositionFactory.eigSymm(N, false)
      decomposer.decompose(denseL)

      // obtain the eigenvalues and compute a histogram of the results
      val devals = (0 until N).map{ i => decomposer.getEigenvalue(i).real }.toArray.sorted
      val hist = new Histogram(devals, Option(0.0, 2.0), 10)
      val scaleFact = 1.0 / hist.hist.sum // 1.0 / ((2.0 / 100) * hist.hist.sum)
      val hvals = hist.hist.map{ h => h.toDouble * scaleFact }
      //println( hvals.mkString(", ") )
      hvals
    } else {
      return Array.ofDim[Double](0)
    }
  }

}
