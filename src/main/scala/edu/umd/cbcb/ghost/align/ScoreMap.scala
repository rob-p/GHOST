package edu.umd.cbcb.ghost.align

import scala.collection.parallel.{ ParIterable, ParMap => PMap }

case class ScoreMap(val map: Option[PMap[(String,String), Double]], val mv: Double) { 
  private val _map = map.getOrElse( PMap.empty[(String,String), Double] )
  private val _getAndScale = _map
  def apply(k: (String,String)) = { 
    _map.get(k) match {
      case Some(x) => x 
      case None => mv
    }
  }

}
