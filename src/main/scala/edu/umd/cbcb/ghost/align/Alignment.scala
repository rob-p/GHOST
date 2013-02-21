package edu.umd.cbcb.ghost.align

import scala.collection.JavaConversions._
import scala.collection.mutable.{ OpenHashMap => MMap }
import com.google.common.collect.{ HashBiMap => MBiMap }

object Alignment {

  def empty[N] = { new Alignment[N]() }

  case class AlignmentIterator[N]( a:Alignment[N] ) extends Iterator[(N,N,Double)] { 
    private val _alnIt = a._aln.entrySet().iterator()

    def hasNext = { _alnIt.hasNext() }

    def next = {
      val kv = _alnIt.next()
      val (k,v) = (kv.getKey(), kv.getValue())
      val s = a._scores((k,v))
      (k,v,s)
    }
  }

}

class Alignment[N] extends Iterable[(N,N,Double)] {
  // To access the companion's methods
  import Alignment._

  // The bi-directional map that will hold the alignment
  private val _aln:MBiMap[N,N] = MBiMap.create()
  // The map that holds the score for each aligned pair
  private val _scores = MMap.empty[(N,N),Double]

  // Our iterator 
  override def iterator = { AlignmentIterator(this) }

  // The size of the alignment
  override def size() = { _aln.size() }

  def alignedLeft() = { _aln.keySet() }
  def alignedRight() = { _aln.values() }

  // Is lk (from the left set) currently aligned
  def isAlignedLeft( lk: N ) = { _aln.containsKey(lk) }
  // Is rk (from the right set) currently aligned
  def isAlignedRight( rk: N ) = { _aln.containsValue(rk) }
  // Is either lk or rk currently aligned
  def isEitherAligned( lk:N, rk:N ) = { (isAlignedLeft(lk) || isAlignedRight(rk)) }

  def ->-( lk:N ): Option[N] = { 
    if ( isAlignedLeft(lk) ) { 
      Some(_aln.get(lk))
    } else { 
      Option.empty[N]
    }
  }

  def getScore( lk:N, rk:N ) = { 
    _scores.getOrElse( (lk,rk), 0.0 )
  }

  /** Unsafe forward lookup
  *
  * Returns the mapping, in the forward direction, for lk.
  * If no mapping exists; it throws an exception.
  * 
  */
  def ->-!( lk:N ): N = { _aln.get(lk) }
  def -<-!( rk:N ): N = { _aln.inverse().get(rk) }

  def -<-( rk:N ): Option[N] = { 
    if ( isAlignedRight(rk) ) { 
      Some(_aln.inverse().get(rk))
    } else { 
      Option.empty[N]
    }
  }

  def remove( lk: N, rk: N ) = { 
    _scores -= ((lk,rk))
    _aln.remove(lk)
  }

  /** Forcibly inserts the mappking (lk->rk) into the alignment
  *
  * This method will insert (lk->rk) into the alignment even if it
  * causes other nodes to become unaligned.
  *
  * @return An option containing the removed left key; None if no left key
  * was removed.
  */
  def forceInsert( lk: N, rk: N, s: Double ) = { 

    ->-(lk) match { 
      case Some(rkp) => _scores -= ((lk, rkp))
      case None => ()
    }

    val mappedKey = -<-(rk) match { 
      case Some(lkp) => { _scores -= ((lkp,rk)); Some(lkp) }
      case None => { Option.empty[N] }
    }

    _scores((lk,rk)) = s
    _aln.forcePut(lk, rk) 
    mappedKey
  }

  /** Insert the match (lk,rk,s) into the alignment
  *
  * Assumes that neither lk nor rk are currently aligned
  */
  def insert( lk: N, rk: N, s: Double ) = {
    if (!isEitherAligned(lk,rk)){
      _aln.put(lk,rk)
      _scores((lk,rk)) = s
      true
    } else { 
      false
    }
  }

  /** Insert the match (lk,rk,s) into the alignment
  *
  * Convenience method so that the caller doesn't have to
  * unwrap the match. Assumes that neither lk nor rk are
  * currently aligned
  */
  def insertMatch( m: (N,N,Double) ) = {
    val (lk,rk,s) = m
    insert(lk,rk,s)
  }

  def insertIfBetter( lk:N, rk:N, s:Double ) = {
    val matchCheck = (isAlignedLeft(lk), isAlignedRight(rk))

    /*
    * See if the proposed match is better than the
    * current matches.
    */
    matchCheck match {
      case (true, true) => {
        val (rkp, lkp) = (_aln.get(lk), _aln.inverse().get(rk))
        if ((s < _scores((lk,rkp))) && (s < _scores((lkp,rk)))) {
          _aln.forcePut(lk,rk)
          _scores -= ((lk,rkp))
          _scores -= ((lkp,rk))
          _scores += ((lk,rk)) -> s
          true
        } else {
          false
        }
      }
      case (true, false) => {
        val rkp = _aln.get(lk)
        if (s < _scores((lk,rkp))) {
          _aln.forcePut(lk,rk)
          _scores -= ((lk,rkp))
          _scores += ((lk,rk)) -> s
          true
        } else {
          false
        }
      }
      case (false, true) => {
        val lkp = _aln.inverse().get(rk)
        if (s < _scores((lkp,rk))) {
          _aln.forcePut(lk,rk)
          _scores -= ((lkp,rk))
          _scores += ((lk,rk)) -> s
          true
        } else {
          false
        }
      }
      case (false, false) => { 
        _aln.put(lk,rk)
        _scores += ((lk,rk)) -> s
        true
      }
    }

  }
  
}
