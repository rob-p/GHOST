package edu.umd.cbcb.ghost.align.matchers

import scala.math._
import scala.collection._
import scala.collection.parallel.mutable.ParArray
import scala.collection.parallel.ParMap
import scala.collection.parallel.ParIterable

import edu.umd.cbcb.ghost.align.SignatureMap
import org.ejml.simple.SimpleMatrix

class LAPMatcher ( matches: collection.Map[String, ParArray[(String, String, Double)]], ls: SignatureMap, rs: SignatureMap ) extends Matcher(matches, ls, rs) { 

  override def performMatching: ParIterable[MatchType] = { 
    val x = matches.flatMap{ x => x._2 }.map{ kv => 1.0 - kv._3 }.toArray.par
    val xstar = binarizeSolutionOptimal(x, matches.flatMap{ x => x._2 }.toArray.par )
    val retMatch = matches.flatMap{ x => x._2 }.zipWithIndex.filter{ mi => xstar(mi._2) == 1}.map{ mi => mi._1 }
    retMatch.toArray.par
  }

	
}
