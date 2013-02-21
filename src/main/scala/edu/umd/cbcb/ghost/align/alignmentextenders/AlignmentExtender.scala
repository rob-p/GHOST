package edu.umd.cbcb.ghost.align.alignmentextenders

import edu.umd.cbcb.ghost.align.Alignment
import akka.actor.ActorSystem
import net.robpatro.utils.console.ProgressBar
import collection.{ Map => CMap, LinearSeq }
import scala.collection.JavaConversions._
import scala.collection.parallel.mutable.ParArray
import scala.collection.parallel.ParIterable
import collection.mutable.{ Map => MMap, Set => MSet, ConcurrentMap => CMMap, HashMap => MHashMap, HashSet => MHashSet }
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
import org.ejml.simple._
import org.ejml.data._
import org.ejml.ops.NormOps

import collection.parallel.mutable.{ ParSet => MParSet }
import collection.mutable.{ Map => MMap }
import edu.umd.cbcb.ghost.graphtools.Utils
import edu.umd.cbcb.ghost.align.matchers._
import edu.umd.cbcb.ghost.mathutils.TFIDFWeighter
import edu.umd.cbcb.ghost.mathutils.UtilFunctions

class AlignmentExtender(
  val seed: (String, String, Double),
  val ls: SignatureMap,
  val rs: SignatureMap,
  val blastScores: Option[ScoreMap],
  val a: Double,
  val matchType: String,
  val k: Int,
  val ematches: scala.collection.Map[(String, String),Double],
  val alignment: Alignment[SimpleVertex],
  val pbar: ProgressBar,
  val matcherCreator: Option[MatcherCreator],
  val system: ActorSystem ) extends MatchTypes {

  implicit val sys = system
 
  def extendAlignment: collection.Set[(String,String,Double)] = {
    // elements of the queue are ordered by score
    val ord = Ordering[Double].on[MatchType]( m => 1.0 - m._3 )
    val pq = new PriorityQueue[MatchType]()(ord)
    // val pq = MParSet.empty[ MatchType ]
    pq += seed
    var couldMap = true
    val nhop = 3
    val nbins = ls.binWidths.slice(0,nhop).sum

    /* Weight Vector */
    /*
    val avgDeg = 0.5 * (ls.avgDeg + rs.avgDeg)
    var invFeatWeight = 0.0
    val wvecDat = ls.binWidths.zipWithIndex.flatMap{
      wi =>
	val (width, index) = wi
      val fw = 1.0 // pow( avgDeg, -index)
      invFeatWeight += fw
      Array.fill(width){ fw }
    }.toArray
    invFeatWeight = (0.5 / invFeatWeight)
    val wvec = new SimpleMatrix( Array(wvecDat) )
    // var cmatches = ematches // + ((seed._1, seed._2))
    // var mappedLeft = matchedLeft
    // var mappedRight = matchedRight
    // cmatches.map{ x => x._1 } // (x._1, (x._2, x._3)) }
    // cmatches.map{ x => x._2 } // (x._2, (x._1, x._3)) }
    */

    val haveBlastScores = blastScores.isDefined
    var alpha = if(haveBlastScores) { a } else { 1.0 }
    val seqDist = blastScores.getOrElse( ScoreMap(None, 0.0) )

    val vertsG = ls.G.vertexSet
    val vertsH = rs.G.vertexSet

    val avgG = ls.avgDeg
    val avgH = rs.avgDeg

    /*
    val A = new SimpleMatrix( nbins, nbins )

    var offset = 0
    val r = 4.0
    ls.binWidths.slice(0, nhop).map{ w =>
      for( i <- offset until (offset+w);
	   j <- offset until (offset+w) ) {

	     val fw = pow( avgDeg, -w)
	     val d = fw * min( abs(i-j), r )
	     A.set(i, j, 1.0 - (abs(i-j) / w.toFloat))
      }
      offset += w
    }
    */

    val maxGap = 1

    while ( ! pq.isEmpty ) {
      // Take the highest scoring pair from the priority queue
      val (cl1, cr1, s1) = pq.dequeue
      val (clv, crv) = (ls.nameVertexMap(cl1), rs.nameVertexMap(cr1))
      //println("local seed ("+cl1+", "+cr1+"); inserted = "+extendFromSeed)

      /*
      //val a = ((mappedLeft contains cl1) || (mappedRight contains cr1))
      val b = alignment.isEitherAligned( ls.nameVertexMap(cl1), rs.nameVertexMap(cr1))
      //assert(a==b)
      if ( ! b ) { 
	//println("Mapping %s to %s and examining neighborhood".format(cl1, cr1))
        //cmatches += ((cl1, cr1, s1))
        //val currentScore = currentMap.getOrElse(ls.nameVertexMap(cl1), (cr1, Double.PositiveInfinity))._2
        //if ( currentScore > s1 ) { 
        alignment.insert(clv, crv, s1)
	//currentMap( ls.nameVertexMap(cl1) ) = (rs.nameVertexMap(cr1), s1)
        //}
	//mappedLeft += cl1
	//mappedRight += cr1
      }
      */

      alignment.insertIfBetter(clv, crv, s1)
      val extendFromSeed = true
      pbar.update(alignment.size)

      // Get the neighborhood of the paired proteins in their respective networks
      var adjLeft = ls.neighborhoodByName(cl1, 1)//.filterNot{ m => mappedLeft contains m }
      var adjRight = rs.neighborhoodByName(cr1, 1)//.filterNot{ m => mappedRight contains m }

      // which side of the alignment will we consider expanding to account for gaps
      object GapChoice extends Enumeration{
	val Left, Right, Neither = Value
      }

      val expandSide = (adjLeft.size - adjRight.size) match {
	case x: Int if x > 0 => GapChoice.Right // the left is bigger than the right, introduce gaps on the right
	case x: Int if x < 0 => GapChoice.Left // the right is bigger than the left, introduce gaps on the left
	case _ => GapChoice.Neither // they are the same size; introduce no gaps
      }
      
      var allMapped = false
      var consideredGapSize = 0

      while ( extendFromSeed && !allMapped && consideredGapSize < maxGap ) {

	// Compute the pairwise scores
	def computePairwiseScores( adjLeft: collection.Set[SimpleVertex], adjRight: collection.Set[SimpleVertex] ) = {
	  //val idf = TFIDFWeighter.computeIDFVector( ls.dmat, rs.dmat )

	  Random.setSeed(System.currentTimeMillis)

	  val matches = adjLeft.par.flatMap{
	    lv =>
	      val lname = lv.name
	      //val lsig = ls.dmat.extractVector(true, ls.nameToId(lname))
	      //val ldeg = ls.G.degreeOf( ls.nameVertexMap(lname) )

	      val ln = ls.neighborhoodByName(lname, 1).toSet

              //var lvec = lsig.elementMult(idf)
              var lspec = ls.spectra(lname)
              // lvec = lsig.scale( 1.0 / NormOps.normP2( lsig.getMatrix ) )

	      val dists = adjRight.collect{
                //case rv if( seqDist(lname, rv.name) <= 1.0 ) =>
                case rv =>
	          val rname = rv.name                

                  if ( ematches contains (lname,rname) ) { 
                    
                    val dist = ematches((lname, rname))
                    if ( alignment->-(lv) == rv ){ 
                      (rname, 0.0)
                    } else { 
                      (rname, dist)
                    }
                  } else { 

	            //val rn = rs.neighborhoodByName(rname,1).toSet
	            //val rdeg = rs.G.degreeOf( rs.nameVertexMap(rname) )
	            //val rsig = rs.dmat.extractVector(true, rs.nameToId(rname))//.extractMatrix(0,1,0,nbins)
                    //var rvec = rsig.elementMult(idf)

                    var rspec = rs.spectra(rname)

                    val structDist = UtilFunctions.structDistance( lspec, rspec )
                    assert(structDist < 1.0)
                    var dist = alpha * structDist

	            //var dist = alpha * (invFeatWeight * UtilFunctions.l1Distance( wvec, lvec, rvec ) )
                    // var dist = alpha * UtilFunctions.l1Distance( wvec, lvec, rvec )
  	            // var dist = alpha * UtilFunctions.cosineDistance( lvec, rvec )
	            // var dist = alpha * NormOps.normP1( lvec.minus(rvec).getMatrix )
	            // var dist = alpha * NormOps.normP1( lvec.minus(rvec).getMatrix )

	            dist += (1.0 - alpha) * seqDist( (lname, rname) )

	            (rname, dist)
                 }
	    }.toArray

	    val ndists = dists.sortBy{ x => x._2 }.slice(0,k).map{ m => (lname, m._1, m._2) }
	    ndists.par
	  }.seq

	  Random.shuffle(matches)
	  matches
	}

	val pmatches = computePairwiseScores( adjLeft, adjRight ).seq

	// don't consider things already mapped
	val newmatches = pmatches.filterNot{ case (lk,rk,d) => 
          //alignment.isEitherAligned( ls.nameVertexMap(lk), rs.nameVertexMap(rk) )
          //}
	  //m =>
          val b = alignment.isEitherAligned(ls.nameVertexMap(lk), rs.nameVertexMap(rk) )
	  //val a = ((mappedLeft contains lk)  || (mappedRight contains rk))
          //assert(b == a)
          b
	}

	var assignedMatches : Iterable[(String, String, Double)] = Iterable.empty

        if (matcherCreator.isDefined) {
	  assignedMatches = if ( newmatches.size > 1 ) {
	    val lm = matcherCreator.get.create(newmatches.groupBy{ m => m._1 }.map{ x => (x._1, x._2.toArray.par) }, ls, rs)
	    val amatches = lm.performMatching
	    amatches.seq
	  } else {
	    newmatches
	  }
          
          assignedMatches = assignedMatches.filter{ case (lm, rm, d) => seqDist((lm,rm)) <= 1.0 }

	  assignedMatches.foreach{ nm =>
            //case (lk, rk, d) =>
            //val (lkv, rkv) = (ls.nameVertexMap(lk), rs.nameVertexMap(rk))
           
            //cmatches += nm
            val (cl1, cr1, s1) = nm
            //val currentScore = currentMap.getOrElse(ls.nameVertexMap(cl1), (cr1,Double.PositiveInfinity))._2
            //if ( currentScore > s1 ) { 
            alignment.insert(ls.nameVertexMap(nm._1), rs.nameVertexMap(nm._2), nm._3)
	    //currentMap( ls.nameVertexMap(cl1) ) = (rs.nameVertexMap(cr1), s1)
            //}
	    //currentMap( ls.nameVertexMap(nm._1) ) = (rs.nameVertexMap(nm._2), nm._3)
	    //mappedLeft += nm._1
	    //mappedRight += nm._2

	  }
	} else {
	  assignedMatches = newmatches
	}

	// expand the appropriate side
	consideredGapSize += 1
	expandSide match {

	  case GapChoice.Left => {
	    adjRight = adjRight.filterNot{ v => alignment.isAlignedRight(v) }//mappedRight contains v }
	    adjLeft = (adjLeft ++ ls.neighborhoodByName(cl1, consideredGapSize+1)).filterNot{ m => alignment.isAlignedLeft(m) }//mappedLeft contains m }
	    allMapped = (adjRight.size == 0)
	  }

	  case GapChoice.Right => {
	    adjLeft = adjLeft.filterNot{ v => alignment.isAlignedLeft(v) } //mappedLeft contains v }
	    adjRight = (adjRight ++ rs.neighborhoodByName(cr1, consideredGapSize+1)).filterNot{ m => alignment.isAlignedRight(m) }//mappedRight contains m }
	    allMapped = (adjLeft.size == 0)
	  }

	  case GapChoice.Neither => {
	    allMapped = true
	  }

	}

	// Place the new potential seeds on the priority queue
	if (consideredGapSize <= 1) {
	  pq ++= assignedMatches.map{ t => (t._1, t._2, t._3)}
	}

      }
    }

    Set.empty[(String, String, Double)]
  }

}


/*
  def retrievePairwiseScores( adjLeft: collection.Set[SimpleVertex], adjRight: collection.Set[SimpleVertex] ) = {
  	  Random.setSeed(System.currentTimeMillis)
  	  val matches = adjLeft.par.flatMap{ lv =>
  	    val lname = lv.name

  	    val dists = adjRight.map{ rv =>
  	      val rname = rv.name
  	      ( lname, rname, ematches( (lname, rname) ) )
  	    }
  	    dists
    	  }.seq
  	matches
      }

*/

// excessLeft.foreach{ x => x.updateAge }; excessRight.foreach{ x => x.updateAge }
// val (remLeft, remRight) = (excessLeft.filter{ x => (x.age >= maxGap) || (mappedLeft contains x.v) },
// 				 excessRight.filter{ x => (x.age >= maxGap) || (mappedRight contains x.v) })

// 	adjLeft.filterNot{ v => (mappedLeft contains v) || (remLeft contains v) }.foreach{ v => excessLeft += new ExcessVertex(v) }
// 	adjRight.filterNot{ v => (mappedRight contains v) || (remRight contains v) }.foreach{ v => excessRight += new ExcessVertex(v) }

// excessLeft = excessLeft &~ remLeft; excessRight = excessRight &~ remRight
