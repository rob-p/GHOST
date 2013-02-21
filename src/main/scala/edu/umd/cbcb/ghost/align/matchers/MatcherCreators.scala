package edu.umd.cbcb.ghost.align.matchers

import akka.actor.ActorSystem
import scala.collection.parallel.mutable.ParArray
import scala.collection.parallel.immutable.ParMap
import scala.collection.parallel.ParIterable
import collection.mutable.{ ConcurrentMap => CMMap, HashMap => MHashMap, HashSet => MHashSet, PriorityQueue }
import edu.umd.cbcb.ghost.align.SignatureMap

/** Classes to create different kinds of matchers */
abstract class MatcherCreator {
  def create( matches: collection.Map[String, ParArray[(String, String, Double)]], ls: SignatureMap, rs: SignatureMap ) (implicit system:ActorSystem) : Matcher
}

object LAPMatcherCreator extends MatcherCreator{
  override def create( matches: collection.Map[String, ParArray[(String, String, Double)]], ls: SignatureMap, rs: SignatureMap )(implicit system:ActorSystem) : Matcher = {
    new LAPMatcher( matches: collection.Map[String, ParArray[(String, String, Double)]], ls: SignatureMap, rs: SignatureMap )
  }
}

object SpectralMatcherCreator extends MatcherCreator{
  override def create( matches: collection.Map[String, ParArray[(String, String, Double)]], ls: SignatureMap, rs: SignatureMap )(implicit system:ActorSystem) : Matcher = {
    new SpectralMatcher( matches: collection.Map[String, ParArray[(String, String, Double)]], ls: SignatureMap, rs: SignatureMap, system: ActorSystem )
  }
}
