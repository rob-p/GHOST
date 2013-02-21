package edu.umd.cbcb.ghost.align.localsearch

import net.robpatro.utils.console.ProgressBar
import edu.umd.cbcb.ghost.io.{ SimpleVertex }
import scala.util.control.Breaks.{ break, breakable }
import edu.umd.cbcb.ghost.align.{ Alignment, ScoreMap }
import scala.collection.JavaConversions._
import edu.umd.cbcb.ghost.graphtools._
import org.jgrapht.graph._
import scala.collection.mutable.{ HashSet => MSet }

class NaiveSearch[V <: SimpleVertex, E, G[V,E] <: Pseudograph[V,E]]( 
  aln: Alignment[V], 
  G : G[V,E], 
  H : G[V,E], 
  blastScores: Option[ ScoreMap ] ) { 


    type ScoredMatch = (V,V,Double)
    val seqDist = blastScores.getOrElse( ScoreMap(None, 0.0) )

   def searchIter( ncRatio: Double ) = {  //constrained: Boolean = false ) = { 

    //val f = Utils.inducedConservedStructure(aln.->-!, G, H)
    val alignmentSize = aln.size
    var toOptimize = aln.map{ case m => m._1 }.toSet
    
    var i = 0
    var ni = 0

    // The edges of H mapped under the current alignment
    var mappedEdges = Utils.mappedEdges( aln.->-!, G, H )
    var f = mappedEdges.size

    // The mapped vertices under the alignment, and the induced
    // subgraph of H
    println("order(G)="+G.vertexSet.size+", |alignment|="+aln.size)
    val mappedVerts = G.vertexSet.map{ v => aln.->-!(v) }
    val inducedH = new Subgraph[V,E,G[V,E]](H, mappedVerts)

    def isMapped( fm: (V) => V, l: V, r: V ): Int = { 
      if (G.containsEdge(l,r) && H.containsEdge(fm(l), fm(r))) { 1 } else { 0 }
    }


    case class LocalMove( inserts: Seq[ScoredMatch], delta: Int ) { 
      
      def perform() = { 
        var replacedKeys = MSet.empty[V]
        var removedKeys = inserts.flatMap{ case(l,r,s) => 
          replacedKeys += l
          aln.forceInsert(l,r,s) 
        }.toSet
        removedKeys &~ replacedKeys
      }

    }
    
    // Ordering on the move priority queue
    val pqord = Ordering[Double].on[LocalMove]{ m => m.delta }

    val pb = new ProgressBar( toOptimize.size, "@" )
    
    while (!toOptimize.isEmpty) { 
      val l = toOptimize.head
      toOptimize -= l

      // Check if l is currently aligned
      val (r, v) = aln.->-(l) match { 
        // If not, assign it to a random unassigned node in H
        case None => { 
          val tmpr = (H.vertexSet &~ aln.alignedRight()).head
          val score = seqDist((l.name, tmpr.name))
          // insert this temporary match
          aln.insert(l, tmpr, score)
          (tmpr,score)
        }
        // Otherwise, return it's currently aligned node
        case Some(r) => (r, seqDist((l.name, r.name)))
      }
      
      // For each pair currently in the alignment
      // alignmentEdges.foreach{ case (l,r,v) if ( aln.->-!(l) == r ) =>
      // breakable{ 
      // for every nl \in N(l), count the number of 
      // how many times is f(l) = r, f(nl) an edge in H
      val nL = Utils.neighbors(l, G)
      val mappedLR = nL.collect{ case x if (H.containsEdge(aln->-!(x),r)) => Set(aln->-!(x),r) }
      
      var bestMove = LocalMove( Seq.empty[ScoredMatch], 0  )

        val constrained = (scala.util.Random.nextDouble > ncRatio)
      // Check for swapping (l->r) with (l->rp)
      val moves = H.vertexSet.par.map{ rp =>

        var delta = 0
        aln-<-(rp) match { 
         // If rp is not currently mapped
         // (l->r) --> (l->rp)
          case None => { 
            val removed = mappedLR
            val mappedLRp = nL.collect{ case x if (H.containsEdge(aln->-!(x),rp)) => Set(aln->-!(x),rp) }
            val added = mappedLRp
            val deltaEdges = added &~ removed
            delta += deltaEdges.size - (removed &~ added).size

            val sd = seqDist( (l.name, rp.name) )

            var score = if ( delta > bestMove.delta && sd <= v ) { delta } else { -1 }

              //aln.forceInsert(l,rp,0.0)
              //assert( aln.size == alignmentSize, { println("swap "+(l->r)+" --> "+(l->rp)+" changed alignment size from "+
              //alignmentSize+" to "+aln.size)})
            //bestMove = LocalMove( List( (l, rp, sd) ), delta )
            LocalMove( List( (l, rp, sd) ), score )
          }

         // If rp is mapped (lp->rp), perform the swap 
         // (l->r) --> (l->rp) 
         // (lp->rp) --> (lp->r)
          case Some(lp) => { 
            var score = -1
            var moveSet = List.empty[ScoredMatch]
            if (lp != l ) { 
              val nLp = Utils.neighbors(lp, G)
              val mappedLpRp = nLp.collect{ case x if (H.containsEdge(aln->-!(x),rp)) => Set(aln->-!(x),rp) }
              val removed = mappedLR | mappedLpRp
            
              val mappedLRp = nL.collect{ case x if (H.containsEdge(aln->-!(x), rp)) => Set(aln->-!(x),rp) }
              val mappedLpR = nLp.collect{ case x if (H.containsEdge(aln->-!(x), r)) => Set(aln->-!(x),r) }
              val added = mappedLRp | mappedLpR

              val deltaEdges = added &~ removed
              delta += (deltaEdges.size) - (removed &~ added).size
              val sdLRp = seqDist( (l.name,rp.name) )
              val sdLR = seqDist( (l.name,r.name) )
              val sdLpR = seqDist( (lp.name,r.name) )
              val sdLpRp = seqDist( (lp.name,rp.name) )

              val constrainedImprovement = if (constrained) { sdLpR <= sdLpRp } else { true }
              if ( delta > bestMove.delta && sdLRp <= sdLR && constrainedImprovement ) { 
                moveSet = List( (l, rp, sdLRp), (lp, r, sdLRp) )
                score = delta 
              } else {
                score = -1 
              }
             
                   //seqDist( (lp.name, r.name) ) <= 1.0 ) { 
                //assert( aln->-!(lp) == rp , { println("lp -/-> rp") })
                //assert( aln->-!(l) == r ,{ println("l -/-> r") })
                //assert( aln-<-!(rp) == lp ,{ println("lp <-/- rp") })
                //assert( aln-<-!(r) == l ,{ println("l <-/- r") })
                //aln.remove(lp,rp)
                //aln.remove(l,r)
                //aln.insert(lp,r,0.0)
                //aln.insert(l,rp,0.0)
                     //bestMove = LocalMove( List( (lp,r,sdLpR), (l,rp,sdLRp) ), delta )
                //bestMove = LocalMove( List( (l, rp, sdLRp), (lp, r, sdLRp) ), delta )
                //assert( aln.size == alignmentSize, { println("swap "+(l->r)+" --> "+(l->rp)+" and "+
                //(lp->rp)+" --> "+(lp->r)+" changed alignment size from "+alignmentSize+" to "+aln.size)})
            }
                LocalMove( moveSet, score )
            }
          }
        }
       
        
        bestMove = moves.maxBy{ 
          case LocalMove( inserts, delta ) => if ( !inserts.isEmpty ) { (delta, 1.0 - inserts.head._3) } else { (-1, -1.0) }
        }
       /*
       if (delta > 0) { 
         println("improved score by "+delta)
         break()
       }
       */

        //}

       pb.update(i)

       if ( bestMove.delta > 0 ) { 
         val preSize = toOptimize.size
         val removedNodes = bestMove.perform
         removedNodes.foreach{ lp => 
           println("added "+lp.name+" to toOptimize ")
           toOptimize += lp 
         }
         i += 1-(removedNodes.size)
       } else { 
         i += 1
       }
                        //}
                        //}
      
    }
    pb.done()

  }
   
   
  /** Override an alignment with a single pair
  *
  */
  def omap( m: scala.collection.Map[V,V] ) = { 
    z:V => if (m contains z) { m(z) } else { aln.->-!(z) }
  }

}

/*
class MoveEvaluator[V, E, G[V,E] <: Pseudograph[V,E]](
  val G: G[V,E],
  val H: G[V,E],
  ) extends Actor { 

  /*
   * The graph and k hop neighborhood with which
   * this actor will work
   */
  private var constrained = false

  def receive = {
    // Receive a request to evaluate a move
    case EvaluateMoveRequest(fwMap, move) => rec ! evaluateMove(fwMap, move)
    case SetConstrained(cval) => constrained = cval
    case ExtractSubgraphDescriptor(vert) => rec ! SubgraphDescriptorMsg(generateSubgraphDescriptor(vert))
    case d: DeadLetter => println(d)
    case _ => println("received unknown message")
  }

  def evaluateMove( fwMap: (V) => V, proposedMove: (V,V) ) = { 
      val scoredMove = aln-<-(rp) match { 
         // If rp is not currently mapped
         // (l->r) --> (l->rp)
          case None => { 
            val removed = mappedLR
            val mappedLRp = nL.collect{ case x if (H.containsEdge(aln->-!(x),rp)) => Set(aln->-!(x),rp) }
            val added = mappedLRp
            val deltaEdges = added &~ removed
            delta += deltaEdges.size - (removed &~ added).size

            val sd = seqDist( (l.name, rp.name) )

            if ( delta > bestMove.delta && sd <= v ) { 
                   bestMove = LocalMove( List( (l, rp, sd) ), delta )
            }

          }

         // If rp is mapped (lp->rp), perform the swap 
         // (l->r) --> (l->rp) 
         // (lp->rp) --> (lp->r)
          case Some(lp) => { 
            if (lp != l ) { 
              val nLp = Utils.neighbors(lp, G)
              val mappedLpRp = nLp.collect{ case x if (H.containsEdge(aln->-!(x),rp)) => Set(aln->-!(x),rp) }
              val removed = mappedLR | mappedLpRp
            
              val mappedLRp = nL.collect{ case x if (H.containsEdge(aln->-!(x), rp)) => Set(aln->-!(x),rp) }
              val mappedLpR = nLp.collect{ case x if (H.containsEdge(aln->-!(x), r)) => Set(aln->-!(x),r) }
              val added = mappedLRp | mappedLpR

              val deltaEdges = added &~ removed
              delta += (deltaEdges.size) - (removed &~ added).size
              val sdLRp = seqDist( (l.name,rp.name) )
              val sdLR = seqDist( (l.name,r.name) )
              val sdLpR = seqDist( (lp.name,r.name) )
              val sdLpRp = seqDist( (lp.name,rp.name) )

              val constrainedImprovement = if (constrained) { sdLpR <= 1.0 } else { true }
              if ( delta > bestMove.delta && 
                   sdLRp <= sdLR &&
                   constrainedImprovement ) { 
                bestMove = LocalMove( List( (l, rp, sdLRp), (lp, r, sdLRp) ), delta )
              }
            }
          }
        }
  
  }

}
*/
