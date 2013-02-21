package edu.umd.cbcb.ghost.align.matchers

import scala.concurrent.duration._
import scala.concurrent.{ Await, Future, ExecutionContext, Promise }

import akka.actor.{ Actor, ActorSystem, DeadLetter, Props }
import akka.pattern.ask

import com.typesafe.config.ConfigFactory

import scala.math._
import scala.util.{ Random, Try, Success, Failure }
import scala.collection._
import scala.collection.JavaConversions._
import scala.collection.parallel.mutable.ParArray
import scala.collection.parallel.immutable.ParMap
import scala.collection.parallel.ParIterable
import scala.collection.mutable.ListBuffer
import java.util.concurrent.{ ConcurrentHashMap => PHMap }

import edu.umd.cbcb.ghost.io.SimpleVertex
import edu.umd.cbcb.ghost.align.SignatureMap
import edu.umd.cbcb.ghost.matchers.actors._

import org.jgrapht.alg.FloydWarshallShortestPaths

import org.ejml.simple.SimpleMatrix
import org.ejml.data.DenseMatrix64F
import org.ejml.ops.NormOps

import org.apache.commons.math.linear._

import net.robpatro.utils.time.Timer
import net.robpatro.utils.console.ProgressBar
import edu.umd.cbcb.mathtools.EigenDecomposerARPACK
import edu.umd.cbcb.mathtools.EigenDecompositionResult
import edu.umd.cbcb.mathtools.{ MatVecMult, ApacheMatVecMult, COOMatVecMult }
import edu.umd.cbcb.ghost.mathutils.TFIDFWeighter
import edu.umd.cbcb.ghost.mathutils.UtilFunctions

import no.uib.cipr.matrix.sparse.FlexCompRowMatrix
import no.uib.cipr.matrix.AbstractMatrix
import no.uib.cipr.matrix.DenseVector
import no.uib.cipr.matrix.Vector.Norm
import scala.util.Random
import ExecutionContext.Implicits.global

object PowerIteration {
  def leadingEig(A: AbstractMatrix, n: Int, thresh: Double = 1e-3) = {
    var bk = new DenseVector(Array.fill(A.numRows) { 1.0 })
    var bk1 = new DenseVector(A.numRows)
    var i = 0
    var err = Double.MaxValue
    while (i < n && err > thresh) {
      A.mult(bk, bk1)
      val norm = bk1.norm(Norm.Two)
      bk1.scale(1 / norm)
      bk.scale(-1)
      bk.add(bk1)
      err = bk.norm(Norm.Two)
      bk.set(bk1)
      i += 1
    }
    bk
  }
}

class SpectralMatcher(matches: Map[String, ParArray[(String, String, Double)]], ls: SignatureMap, rs: SignatureMap, sys: ActorSystem) extends Matcher(matches, ls, rs) {
  var first = true

  private def computeJointSpectra(S: SignatureMap, nodes: ParArray[String]): PHMap[(String, String), SimpleMatrix] = {
    val N = nodes.size
    var reqNumSpec = 0 //(N * N) / 2
    val cmap = new PHMap[(String, String), SimpleMatrix]

    // Create the actors that will be responsible for computing the joint spectra
    def createJointSpectrumActors(n: Int, k: Int) = {
      val extractorArray = new Array[JointSpectrumComputer](n)
      (0 until n).foreach {
        i =>
          extractorArray(i) = new JointSpectrumComputer(S, k, cmap)
          extractorArray(i).start
      }
      extractorArray
    }

    val extractors = createJointSpectrumActors(30, 2)
    val numExtractors = extractors.size

    val rgen = new Random()
    val nnodes = nodes.size
    for (i <- 0 until nnodes; j <- i until nnodes) {
      reqNumSpec += 1
      extractors(rgen.nextInt(numExtractors)) ! ComputeJointSpectrum(nodes(i), nodes(j))
    }

    val pb = new ProgressBar(reqNumSpec, "*")

    while (cmap.size < reqNumSpec) {
      pb.update(cmap.size)
      Thread.sleep(1500)
    }

    pb.update(cmap.size)
    pb.done

    print("\033[1A\033[1A\033[2K\r")

    (0 until numExtractors).foreach { i => extractors(i) ! Stop() }
    val specMap = cmap
    assert(specMap.size >= reqNumSpec, "SpecMap size is " + specMap.size + ", but we need " + reqNumSpec + " entries ")
    specMap
  }

  private def fsweight(nu: Set[SimpleVertex], nv: Set[SimpleVertex], navg: Double) = {
    val uiv = (nu & nv).size
    val uuv = (nu | nv).size
    val umv = (nu &~ nv).size
    val vmu = (nv &~ nu).size
    val lambda = max(0.0, navg - (umv + uiv))

    ((2.0 * uiv) / (umv + 2.0 * uiv + lambda)) * ((2.0 * uiv) / (vmu + 2.0 * uiv + lambda))
  }

  private def jaccardDistance(nu: Set[String], nv: Set[String]) = {
    val uiv = (nu & nv).size.toFloat
    val uuv = (nu | nv).size.toFloat
    1.0 - (uiv / uuv)
  }

  private def constructQAP() = {
    implicit val system = sys
    var weights = matches.flatMap { kv => kv._2 }.map { kv => 1.0 - kv._3 }
    var scale = 1.0 / weights.max
    weights = weights.map { w => w * scale }

    val kN = weights.size
    val maxval = max(1, (kN * (kN + 1) / 2))
    var i = 0
    val (maxWeight, minWeight) = (weights.max, weights.min)
    val mdiag = weights

    Random.setSeed(System.currentTimeMillis)
    val nodesG = Random.shuffle(matches.keys.toList).sortBy { s: String => ls.avgDensity(s) }.toArray.par
    val nodesH = Random.shuffle(matches.flatMap { x => x._2 }.map { x => x._2 }.toList).sortBy { s: String => rs.avgDensity(s) }.toArray.par

    // val jointSpectraG = computeJointSpectra( ls, nodesG )
    // val jointSpectraH = computeJointSpectra( rs, nodesH )

    scale = 1.0 / kN
    var (ma, mi) = (-Double.MaxValue, Double.MaxValue)

    // Map from each key to the index in the flattened map of matches where it's entries start
    val keyIndMap = matches.keysIterator.drop(1).scanLeft((matches.head._1, 0)) { (x, y) => (y, x._2 + matches(x._1).size) }.toMap

    val vertsG = ls.G.vertexSet
    val vertsH = rs.G.vertexSet
    val navg = ((1.0 / vertsG.size) * vertsG.map { v => ls.G.degreeOf(v) }.sum) +
      ((1.0 / vertsH.size) * vertsH.map { v => rs.G.degreeOf(v) }.sum)

    //val idf = TFIDFWeighter.computeIDFVector( ls.dmat, rs.dmat )

    var ctr = 0

    import scala.actors.{ Futures => SFutures }
    import scala.actors.{ Future => SFuture }

    import scala.collection.mutable.ArrayBuffer

    case class MatrixEntries(ri: ArrayBuffer[Int], ci: ArrayBuffer[Int], e: ArrayBuffer[Double])

    /**
     * Compute the entries of the QA matrix relating to matches
     * involving vertex k \in G.
     *
     * @return an Akka Future holding the MatrixEntries object
     */
    def computeEntries(k: String): Future[MatrixEntries] = {

      Future {
        val ri = ArrayBuffer.empty[Int] //new java.util.concurrent.ConcurrentLinkedQueue[Int]
        val ci = ArrayBuffer.empty[Int] //new java.util.concurrent.ConcurrentLinkedQueue[Int]
        val e = ArrayBuffer.empty[Double] //new java.util.concurrent.ConcurrentLinkedQueue[Double]
        val nei = ls.neighborhoodByName(k, 1).map { x => x.name }.toSet


        val matchKeys = matches.keys

        matches(k).seq.view.zipWithIndex.foreach { outerMatchInd =>
          val (vi, vip, p) = outerMatchInd._1
          val (vertI, vertIp) = (ls.nameVertexMap(vi), rs.nameVertexMap(vip))
          val rowInd = keyIndMap(vi) + outerMatchInd._2
          var matchCtr = 0

          // add the combined topo-sequence cost to the diagonal
          ri.add(rowInd); ci.add(rowInd);
          e.add(1.0 + 5.0 * exp(-p))

          // The neighborhood of vip in H
          val neip = rs.neighborhoodByName(vip, 1).map { x => x.name }.toSet

          // The set of matches for all vj s.t. vj \in N(vi)
          val adjKeys = nei.filter { n => (matchKeys contains n) && keyIndMap(vi) > keyIndMap(n) }.seq
          val numKeys = adjKeys.size

          // If there are any other matches compatible with (vi, vip)
          if (numKeys > 0) {

            // Score each compatible match
            adjKeys.foreach { nkey =>

              // Compute the affinity of match (vi,vj) with match (vip, vjp)
              matches(nkey).seq.view.zipWithIndex.foreach {
                case (matchCost, offset) =>
                  matchCtr += 1
                  // extract the variables from the match
                  val (vj, vjp, q) = matchCost
                  val (vertJ, vertJp) = (ls.nameVertexMap(vj), rs.nameVertexMap(vjp))

                  if (ls.G.getEdge(vertI, vertJ) != null && rs.G.getEdge(vertIp, vertJp) != null) {

                    val colInd = keyIndMap(vj) + offset

                    val nej = ls.neighborhoodByName(vj, 1).map { x => x.name }.toSet
                    val nejp = rs.neighborhoodByName(vjp, 1).map { x => x.name }.toSet

                    // Using Jaccard similarity
                    /*
                val unionG = (nei | nej).size
                val unionH = (neip | nejp).size
                val e1 = if ( unionG == 0 ) { 1.0 } else { (nei & nej).size.toDouble / unionG }
                val e2 = if ( unionH == 0 ) { 1.0 } else { (neip & nejp).size.toDouble / unionH }
                */

                    // Using structural distance
                    val (ls1, ls2) = (ls.spectra(vi), ls.spectra(vj))
                    val (rs1, rs2) = (rs.spectra(vip), rs.spectra(vjp))
                    val e1 = UtilFunctions.structDistance(ls1, ls2)
                    val e2 = UtilFunctions.structDistance(rs1, rs2)

                    // Using histograms
                    /*
                val (l1, r1) = ( ls.dmat.extractVector(true, ls.nameToId(vi)), ls.dmat.extractVector(true, ls.nameToId(vj)) )
	        val (l2, r2) = ( rs.dmat.extractVector(true, rs.nameToId(vip)), rs.dmat.extractVector(true, rs.nameToId(vjp)) )
                val e1 = UtilFunctions.cosineDistance(l1, r1)
	        val e2 = UtilFunctions.cosineDistance(l2, r2)
                */

                    //val neighProd = nej.size * nejp.size 
                    //val s = if (neighProd > 0) { 1.0 / neighProd } else { 1.0 }
                    //val s = 1.0 / (nei.size * neip.size)
                    val s = (1.0 - (e1 - e2)) / numKeys

                    /*
               if ( (e1+e2) > 0.0  ) {
                 exp( - (abs(e1-e2) / (e1 + e2) ) ) 
               } else { 0.0 }
                */

                    if (s > 0.0) {
                      ri.add(rowInd); ci.add(colInd);
                      ri.add(colInd); ci.add(rowInd);
                      e.add(s); e.add(s);
                    }

                  }

              }
            }

          }
        }

        val retr = ri
        val retc = ci
        val rete = e
        MatrixEntries(retr, retc, rete)
      }
    }

    /**
     * Actor that aggregates entries computed in parallel into a COO matrix
     *
     * @constructor create a new matrix aggregator given the row, column, and entry
     *              array buffer and the number of messages to receive
     * @param rib array buffer holding row indices
     * @param cib array buffer holding column indices
     * @param eb arrab buffer holding matrix entries
     * @param N the number of messages this actor will have to aggregate
     *          to form the matrix
     */
    class MatrixAggregatorActor(
      val rib: ArrayBuffer[Int],
      val cib: ArrayBuffer[Int],
      val eb: ArrayBuffer[Double],
      val N: Int) extends Actor {

      private var numRowsAdded = 0
      import context._

      def receive = {
        case MatrixEntries(ri, ci, e) => {
          // add the new row & col indices and entries
          rib ++= ri
          cib ++= ci
          eb ++= e
          // update the count of received messages
          numRowsAdded += 1
          // if we're done, then stop the actor
          if (numRowsAdded == N) {
            context.stop(self)
          }
        }
        case d: DeadLetter => println(d)
        case _ => println("received unkown message")
      }
    }

    import scala.collection.mutable.Queue

    val rib = new ArrayBuffer[Int]
    val cib = new ArrayBuffer[Int]
    val eb = new ArrayBuffer[Double]
    val pb = new ProgressBar(matches.size, "%")
    println()
    ctr = 0

    val props = Props().withCreator(new MatrixAggregatorActor(rib, cib, eb, matches.size)).withDispatcher("akka.actor.qr-dispatcher")
    val actor = system.actorOf(props, name = "MatrixAggregatorActor")
    val matRows = matches.keysIterator.foreach { k =>

      val f = computeEntries(k)
      f onComplete {
        case Success(r) => {
          actor ! r
          pb.update(ctr); ctr += 1;
        }
        case Failure(failure) => {
          println("failed to compute QAP matrix")
          System.exit(1)
        }
      }

    }

    while (!actor.isTerminated) {
      Thread.sleep(100)
    }

    // Check that we didn't mess anything up by doing it in parallel
    assert(rib.size == cib.size, "FINAL : Found rowind.size = " + rib.size + " but cind.size = " + cib.size)
    assert(rib.size == eb.size, "FINAL : Found rowind.size = " + rib.size + " but e.size = " + eb.size)

    // Finish the progress bar
    pb.done

    // Go back two lines
    print("\033[1A\033[1A")

    scale = 1.0 / ma

    // Scale the matrix
    //val mit = M.iterator(true, 0, 0, kN-1, kN-1)

    val M = new FlexCompRowMatrix(kN, kN)
    (0 until rib.size).foreach {
      i => M.set(rib(i), cib(i), eb(i))
    }
    //M


    // return the matrix
    var ret : MatVecMult = null
    if ( kN > 20 ) {
      ret = new COOMatVecMult(rib.toArray, cib.toArray, eb.toArray, kN, kN)
    } else {
      val M = new OpenMapRealMatrix(kN, kN)
      (0 until rib.size).foreach{
	i => M.setEntry(rib(i), cib(i), eb(i))
      }
      ret = new ApacheMatVecMult(M)
    }
    
    ret

  }

  override def performMatching: ParArray[MatchType] = {

    val M = constructQAP

    val N = matches.flatMap { kv => kv._2 }.map { kv => 1.0 - kv._3 }.size
    // matches.size // * matches.head._2.size
    var evec = Array.empty[Double]

    // With anything but a tiny matrix, use ARPACK
    if (N > 20) {

      val n = min(8,N-1)
      val eigenDecomp = EigenDecomposerARPACK.decompose(M, n, "LM")
      if ( !eigenDecomp.getSuccess ) {
	println("EigenDecomposition Failed!")
      }

      val evals = eigenDecomp.getEvals
      val numEV = evals.size
      //println("Eigenvalues : "+evals)
      val inds = evals.zipWithIndex.sortBy{ vi => vi._1 }.map{  vi => vi._2 }
      val evecs = eigenDecomp.getEvecs

      require( evals.size == evecs.size, println("Size of evals is %d, but size of evecs is %d".format(evals.size, evecs.size) ) )
      //println( evecs.size() )
      //println( inds.mkString(", "))
      evec = evecs( inds.last )

    } else {
      // Otherwise use a simpler Eigendecomposition routine
      val mvm = M.asInstanceOf[ApacheMatVecMult[OpenMapRealMatrix]]
      //val decomposer = new EigenDecompositionImpl( mvm.getMat, 1e-20 )
      val decomposer = new EigenDecompositionImpl( mvm.getMat.add( mvm.getMat.transpose() ).scalarMultiply( 0.5 )  , 1e-100)
      val evals = decomposer.getRealEigenvalues()
      //println("Eigenvalues : {"+evals.mkString(",")+"}")
      val numEV = evals.size
      val inds = evals.zipWithIndex.sortBy{ vi => vi._1 }.map{  vi => vi._2 }
      evec = decomposer.getV.getColumn( inds.last )

    }


    /*
    val tol = 1e-20
    evec = PowerIteration.leadingEig(M, 200, tol).getData
    */
    var x = evec.toArray.par


    //(0 until N).foreach{ i => x(i) = exp(- abs(evecs(i)) ) }
    //while( mit.hasNext ) { x(i) = exp(- abs(mit.next.getValue) ); i += 1 }
    //val x = evec.iterator(true,0,0,N,0).map{ x => exp(-x) }

    // flatten the set of matches to an array of (String, String, Double)
    var x2 = matches.flatMap { x => x._2 }.toArray
    // binarize the matching with dist(v) = 1.0 - evec(v)

    //val xstar = binarizeSolutionOptimal(  x2.map{ case(lm, rm, s) => 1.0 - s }.toArray.par , x2.toArray.par )

    val xstar = binarizeSolutionOptimal(x, x2.toArray.par)

    // extract the matches in the optimal solution and return them
    val retMatch = x2.view.zipWithIndex.collect { case (matchEntry, index) if (xstar(index) == 1) => matchEntry }

    // return eigenvector scores
    // val retMatch = x2.view.zipWithIndex.collect{ case(matchEntry, index) if (xstar(index) == 1) => { val (lm,rm,s) = x2(index); (lm,rm,1.0-s) } }

    retMatch.toArray.par
  }

}
