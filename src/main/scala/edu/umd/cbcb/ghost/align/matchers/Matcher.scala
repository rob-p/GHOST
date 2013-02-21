package edu.umd.cbcb.ghost.align.matchers
import collection.{ Map => CMap, LinearSeq }
import scala.collection.JavaConversions._
import scala.collection.parallel.mutable.ParArray
import scala.collection.parallel.immutable.ParMap
import scala.collection.parallel.ParIterable
import collection.mutable.{ OpenHashMap => MHashMap, HashSet => MHashSet, PriorityQueue }
import collection.concurrent.{Map => CMMap}
import scala.math._
import scala.util.Random

import java.util.concurrent.{ ConcurrentHashMap => PHMap }

import edu.umd.cbcb.ghost.align.SignatureMap
import org.ejml.simple.SimpleMatrix


abstract class Matcher ( val matches: collection.Map[String, ParArray[(String, String, Double)]], val ls: SignatureMap, val rs: SignatureMap ) extends MatchTypes {

  def performMatching: ParIterable[MatchType]

  def binarizeSolutionGreedy  (x: ParArray[Double], matches: ParArray[MatchType] ) = {
    var mapped = MHashSet.empty[String]
    // var matchesWithInd = matches.seq.zipWithIndex
    var pmatches = matches.seq.zipWithIndex.map{ mi => ((mi._1._1, mi._1._2, x(mi._2)), mi._2) }
    val pqord = Ordering[Double].on[((String,String,Double),Int)]{ m => 1.0 - m._1._3 }
    val pq = new PriorityQueue()(pqord)
    pq ++= pmatches

    val matchVec = Array.fill(matches.size){ 0 }

    while (! pq.isEmpty) {
      val (bm, ind) = pq.dequeue
      val (l,r,s) = bm
      if ( ! ((mapped contains "l"+l) || (mapped contains "r"+r)) ) {
	matchVec(ind) = 1
	//println("adding "+l+" -> "+r+", score = "+s)
	mapped += "l"+l
	mapped += "r"+r
      }
    }

    matchVec
  }


  def binarizeSolutionOptimal (x: ParArray[Double], matches:  ParArray[MatchType] ) = {
    // We solve the Minimum Weight Matching, so invert all weights
    val xs = x.map{ x => 1.0 - x }

    // make it easy to map from a node name to an index in our matching matrix and back
    // val assign every node in the first network a unique ID
    Random.setSeed(System.currentTimeMillis)
    val nodesG = Random.shuffle(matches.map{ s => s._1 }.toSet.toList).sortBy{ s: String => ls.avgDensity(s) }
    val nodesH = Random.shuffle(matches.map{ s => s._2 }.toSet.toList).sortBy{ s: String => rs.avgDensity(s) }

    val nodeToIndG = new MHashMap[String,Int]; nodesG.foreach{ n => if(! nodeToIndG.contains(n)) { nodeToIndG(n) = nodeToIndG.size } }
    val nodeToIndH = new MHashMap[String,Int]; nodesH.foreach{ n => if(! nodeToIndH.contains(n)) { nodeToIndH(n) = nodeToIndH.size } }

    //println(nodeToIndG.size+" :: "+nodeToIndH.size)
    val indToNodeG = nodeToIndG.map{ kv => (kv._2, kv._1) }.toMap
    val indToNodeH = nodeToIndH.map{ kv => (kv._2, kv._1) }.toMap

    val pairIndDict: CMMap[(Int,Int), Int] = new PHMap[(Int,Int), Int]

    val lG = nodeToIndG.size
    val lH = nodeToIndH.size

    /* If the networks are of size N and M respectively, our matching matrix
     *  will be N x M.  However, this implementation of the Hungarian algorithm
     *  requires a square matrix, so we'll give it an N x N matrix assuming N >= M w.l.o.g.
     */
    val N = max(lG,lH)
    val M = min(lG,lH)
    val matXS = Array.fill(N,N){ 100.0 }

    /* To use the C library version */
    /*
    import scala.math._
    import com.sun.jna.Pointer;
    import com.sun.jna.Memory;

    val default = Array.fill(N){ 1000.0 }
    val matXS = Array.ofDim[Pointer](N)
    (0 until N).foreach{ i => matXS(i) = new Memory(N*8); matXS(i).write(0, default, 0, N)  }
    */

    // For each match,
    matches.zipWithIndex.foreach{ mi: (MatchType, Int) =>
      val ((i,j,v),k) = mi
      matXS(nodeToIndG(i))(nodeToIndH(j)) =  xs(k)
      // for C version
      //matXS(nodeToIndG(i)).setDouble(8*nodeToIndH(j), x(k))
      pairIndDict( (nodeToIndG(i), nodeToIndH(j)) ) = k
    }

    val (mrows, mcols) = (Array.ofDim[Int](N), Array.ofDim[Int](N))

    /*
    import kostas.HungarianAlgorithm
    val assignment = HungarianAlgorithm.hgAlgorithm(matXS, "min")
    val matchedPairs = (0 until assignment.size).map{ i => (assignment(i)(0), assignment(i)(1)) }
    */

    /*
    // for C version
    //import asp._
    //ASPLib.INSTANCE.asp(N, matXS, mrows, mcols)
    */

    import javaasp.LinearAssignmentSolver
    val aspi = new LinearAssignmentSolver
    aspi.solve(N, matXS, mrows, mcols)
    val matchedPairs = (0 until N).zip(mrows).par

    val newMatches = matchedPairs.filter{ p => pairIndDict.contains(p) }

    val res = newMatches.map{ m => pairIndDict(m) }
    //println("Found "+res.size+" matches ")
    val matchVec = Array.ofDim[Int]( matches.size )
    for (r <- res) { matchVec(r) = 1 }
    matchVec
  }

}
