package edu.umd.cbcb.ghost.align.alignmentextenders

import collection.{ Map => CMap, LinearSeq }
import scala.collection.JavaConversions._
import scala.collection.parallel.mutable.ParArray
import scala.collection.parallel.ParIterable
import collection.mutable.{ ConcurrentMap => CMMap, HashMap => MHashMap, HashSet => MHashSet }
import scala.math._
import scala.util.Random

import java.util.concurrent.{ ConcurrentHashMap => PHMap }

// we need access to the matchers and signatures
import edu.umd.cbcb.ghost.io.SimpleVertex
import edu.umd.cbcb.ghost.align.matchers._
import edu.umd.cbcb.ghost.align.ScoreMap
import edu.umd.cbcb.ghost.align.SignatureMap
import org.ejml.simple.SimpleMatrix
import org.jgrapht.alg.BellmanFordShortestPath
import scala.collection.mutable.PriorityQueue
import org.ejml.simple.SimpleMatrix
import org.ejml.ops.NormOps

import collection.parallel.mutable.{ ParSet => MParSet }
import collection.mutable.{ Map => MMap }
import edu.umd.cbcb.ghost.graphtools.Utils

class ProximityWalkAlignmentExtender( val seed: (String, String, Double), val ls: SignatureMap, val rs: SignatureMap, val blastScores: Option[ScoreMap], val a: Double, val matchType: String, val k: Int, val ematches: scala.collection.mutable.Set[(String, String, Double)], val currentMap: scala.collection.mutable.Map[SimpleVertex, SimpleVertex]) extends MatchTypes { 

  private def jaccardDistance( nu: Set[SimpleVertex], nv: Set[SimpleVertex] ) = { 
    val uiv = (nu & nv).size.toFloat
    val uuv = (nu | nv).size.toFloat
    if (uuv > 0) { 
      1.0 - ( uiv / uuv )
    } else { 
      1.0
    }
  }

  private def fsweight( nu: Set[SimpleVertex], nv: Set[SimpleVertex], navg: Double ) = { 
    val uiv = (nu & nv).size
    val uuv = (nu | nv).size
    val umv = (nu &~ nv).size
    val vmu = (nv &~ nu).size
    val lambda = max( 0.0, navg - (umv + uiv))

    ( (2.0*uiv) / (umv + 2.0*uiv + lambda) ) * ( (2.0 * uiv) / (vmu + 2.0*uiv + lambda))
  }

  private def computePairwiseScores( adjLeft: Set[SimpleVertex], adjRight: Set[SimpleVertex] ) = { 

  }

  def extendAlignment: collection.Set[(String,String,Double)] = ??? 
  /* currently unimplemented / DEPRECATED 
  { 
    // elements of the queue are ordered by score
    val ord = Ordering[Double].on[MatchType]( m => 1.0 - m._3 )
    val pq = new PriorityQueue[MatchType]()(ord)
    // val pq = MParSet.empty[ MatchType ]
    pq += seed
    var couldMap = true
    val nhop = 4
    val nbins = ls.binWidths.slice(0,nhop).sum
    

    /* Weight Vector */
    val avgDeg = 0.5 * (ls.avgDeg + rs.avgDeg)
    var invFeatWeight = 0.0
    val wvecDat = ls.binWidths.slice(0,nhop).zipWithIndex.flatMap{ 
      wi =>
	val (width, index) = wi
      val fw = pow( avgDeg, -index)
      invFeatWeight += fw
      Array.fill(width){ fw }
    }.toArray
    invFeatWeight = (0.5 / invFeatWeight)
    val wvec = new SimpleMatrix( Array(wvecDat) )
      
    var cmatches = ematches.clone() // + ((seed._1, seed._2))
    
    val preMappedLeft = cmatches.map{ x => x._1 } 
    val preMappedRight = cmatches.map{ x => x._2 }
    val (sizebl, sizebr) = (preMappedLeft.size, preMappedRight.size)

    var mappedLeft = preMappedLeft.clone()
    var mappedRight = preMappedRight.clone()

    val haveBlastScores = blastScores.isDefined 
    val alpha = if(haveBlastScores) { a } else { 1.0 }
    val seqDist = blastScores.getOrElse( ScoreMap(None, 0.0) )

    val vertsG = ls.G.vertexSet
    val vertsH = rs.G.vertexSet
    
    val avgG = ls.avgDeg
    val avgH = rs.avgDeg

    class RunningAverage( var avg: Double ) { 
      var cnt = 1
      def apply( newValue: Double ) { 
	avg = ((cnt * avg) + newValue) / (cnt + 1)
	cnt += 1
      }
    }

    val visitedPairs = MHashMap.empty[(String,String), RunningAverage]
    while ( ! pq.isEmpty ) {

      // Take the highest scoring pair from the priority queue
      val (cl1, cr1, s1) = pq.dequeue

      if (! visitedPairs.contains( (cl1, cr1) ) ) {   // ! ((mappedLeft contains cl1) || (mappedRight contains cr1)) ) { 
	println("Mapping %s to %s and examining neighborhood".format(cl1, cr1))

        cmatches += ((cl1, cr1, s1))
	// currentMap( ls.nameVertexMap(cl1) ) = rs.nameVertexMap(cr1)
	visitedPairs += ((cl1, cr1) -> new RunningAverage(s1))
	mappedLeft += cl1 
	mappedRight += cr1 

	println("CURRENT MAPPED = %d".format(cmatches.size))
	// Get the neighborhood of the paired proteins in their respective networks
	println("Examining pair (%s, %s) with score %f".format(cl1, cr1, s1))
	val adjLeft = ls.neighborhoodByName(cl1, 1)
	val adjRight = rs.neighborhoodByName(cr1, 1)

	// Compute the pairwise scores
	def computePairwiseScores( adjLeft: collection.Set[SimpleVertex], adjRight: collection.Set[SimpleVertex] ) = { 
	  adjLeft.par.flatMap{ 
	    lv =>
	      val lname = lv.name
	    val lsig = ls.dmat.extractVector(true, ls.nameToId(lname)).extractMatrix(0,1,0,nbins)
	    val ln = ls.neighborhoodByName(lname, 1).toSet
	    val ldeg = ls.G.degreeOf( ls.nameVertexMap(lname) )

	    val dists = adjRight.map{ rv =>
	      val rname = rv.name
	      val rn = rs.neighborhoodByName(rname,1).toSet

	      val rdeg = rs.G.degreeOf( rs.nameVertexMap(rname) )
	      val rsig = rs.dmat.extractVector(true, rs.nameToId(rname)).extractMatrix(0,1,0,nbins)
	      var dist = alpha * (invFeatWeight * NormOps.normP1( lsig.minus(rsig).elementMult(wvec).getMatrix )) 
	      dist += (1.0 - alpha) * seqDist( (lname, rname) )
	      // dist *= 1.0 - ( 0.5 * (ldeg + rdeg)  / ( ls.maxDeg + rs.maxDeg ).toFloat )
	      (rname, dist)
	    }.toArray

	    val ndists = dists.sortBy{ x => x._2 }.slice(0,k).map{ m => (lname, m._1, m._2) }
	    ndists.par
	  }.toArray.par
	}

	val pmatches = computePairwiseScores( adjLeft, adjRight ).seq

	// don't consider things already mapped
	val newmatches = pmatches.filterNot{ m => 
	    (mappedLeft contains m._1)  || (mappedRight contains m._2)
	}
	
	val tmatches = newmatches.filter{ x => x._3 < 1.1 * s1 }
	cmatches ++= tmatches

	// Place the new potential seeds on the priority queue
	pq ++= newmatches.map{ t => (t._1, t._2, t._3)}
      } else { 
	println("Incrementing count of (%s,%s)".format(cl1, cr1))
	visitedPairs( (cl1, cr1) )( s1 )  
      }
    }
    
    val consideredMatches = visitedPairs.filterNot{ kv => val x = kv._1; (preMappedLeft contains x._1) || (preMappedRight contains x._2) }
    pq.clear()

    val maxQ = new PriorityQueue[MatchType]()
    maxQ ++= consideredMatches.map{ kv => val (pair, score) = kv; (pair._1, pair._2, score.cnt * (1.0 - score.avg)) }
    val rmatches = collection.mutable.HashSet.empty[(String, String, Double)]

    while( ! maxQ.isEmpty ) { 
      val t = maxQ.dequeue
      if ( ! ((preMappedLeft contains t._1) || (preMappedRight contains t._2)) ) { 
	preMappedLeft += t._1; preMappedRight += t._2
	println(t._3)
	rmatches += t
      }
    }

    // val matcher = new LAPMatcher(consideredMatches.toArray.par, ls, rs)
    // val rmatches = matcher.performMatching
    rmatches.seq.toSet
  }
  */

}
