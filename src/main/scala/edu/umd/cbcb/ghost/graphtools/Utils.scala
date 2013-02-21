package edu.umd.cbcb.ghost.graphtools

import edu.umd.cbcb.ghost.align.{ Alignment }

import scala.collection.JavaConversions._
import scala.math._

import org.ejml.simple._
import edu.umd.cbcb.ghost.io.{ SimpleVertex }
import org.jgrapht.graph._
import org.jgrapht.{ UndirectedGraph }
import org.jgrapht.alg.{ ConnectivityInspector }

import net.robpatro.utils.console.ProgressBar

object Utils {

  def mappedEdges[V, E, G[V,E] <: Pseudograph[V,E]]( fwAln: (V) => V, G: G[V,E], H: G[V,E] ) = { 
    val medge = new PartialFunction[E,E] { 
      def apply( e: E ) = { 
        val (s,t) = (G.getEdgeSource(e), G.getEdgeTarget(e))
        val (sp, tp) = (fwAln(s), fwAln(t))
        H.getEdge(sp,tp)
      }
  
      def isDefinedAt( e: E ) = { 
        val (s,t) = (G.getEdgeSource(e), G.getEdgeTarget(e))
        val (sp, tp) = (fwAln(s), fwAln(t))
        H.containsEdge(sp,tp)
      }
    }

    G.edgeSet.collect{ medge }
  }

  def neighbors[V, E, G[V,E] <: Pseudograph[V,E]]( v: V, G: G[V,E] ) = { 
    G.edgesOf(v).map{ e => 
      val (x,y) = (G.getEdgeSource(e), G.getEdgeTarget(e))
      if (x == v) { y } else { x }
    }
  }

  def inducedConservedStructure[V, E, G[V,E] <: Pseudograph[V,E]]( fwAln: (V) => V, G : G[V,E], H : G[V,E] ) = { 
    
    val mEdges = mappedEdges(fwAln, G, H)
    val mappedVerts = G.vertexSet.collect{ case v if(fwAln(v)!=null) => fwAln(v) }
    val inducedH = new Subgraph[V,E,G[V,E]](H, mappedVerts)
  
    mEdges.size.toDouble / inducedH.edgeSet.size
  }


  def largestConnectedComponent[V, E, G[V,E] <: Pseudograph[V,E]]( g: G[V,E] ): Subgraph[V, E, G[V,E]] = {
    type SubgraphT = Subgraph[V, E, G[V,E]]
    /* Extract the largest connected component from each graph */
    var ci = new ConnectivityInspector(g)
    // Sort the connected components in *descending* order
    var vertexSets = ci.connectedSets.sortWith { (x, y) => x.size > y.size }
    // Grab the largest connected component of G
    if ( vertexSets.size > 0 ) { 
      new SubgraphT(g, vertexSets(0))
    } else { 
      new SubgraphT(g, g.vertexSet)
    }
  }

  def structuralScore[V, E, G[V,E] <: Pseudograph[V,E]]( currentMap: collection.Map[V, V], ng: collection.Set[V], g: G[V,E], h: G[V,E] ) = {
    type SubgraphT = Subgraph[V, E, G[V,E]]
    val mappedNodes = ng.filter{ x => currentMap contains x}.map{ x => currentMap(x) }
    val indSGH = new SubgraphT(h, mappedNodes)
    val indSGG = new SubgraphT(g, ng)

    val indEdgeSize = indSGH.edgeSet.size

    if (indEdgeSize > 0) {
      pow(indSGG.edgeSet.size.toFloat,2) / indSGH.edgeSet.size
    } else {
      0.0
    }
  }


  def otherEndpoint[V <: SimpleVertex, E, G[V,E] <: AbstractBaseGraph[V,E]]( Graph: G[V,E] , c: V, e: E) = {
    val (u,v) = (Graph.getEdgeSource(e), Graph.getEdgeTarget(e))
    val o = if (u == c) { v } else { u }
    o
  }

  def gomoryHuTreeToMatrix[V <: SimpleVertex, E, G[V,E] <: AbstractBaseGraph[V,E]]( ghTree : G[V,E] ) = {
      import org.jgrapht.traverse.BreadthFirstIterator
      // Convert the tree to the max-flow matrix
      val N = ghTree.vertexSet.size
      var n = 0
      val pb = new ProgressBar( N, "*" )
      val maxFlowMat = new SimpleMatrix(N,N)
      maxFlowMat.set(Double.MaxValue)
      ghTree.vertexSet.foreach{
	v =>

	  maxFlowMat.set(v.id, v.id, 100.0 )
	//ghTree.edgesOf(v).foreach{ e => maxFlowMat.set(v.id, otherEndpoint(ghTree,v,e).id, ghTree.getEdgeWeight(e)) }

	  val dfi = new BreadthFirstIterator(ghTree, v)
	  while( dfi.hasNext ) {

	    val nv = dfi.next
	    if ( nv != v ) {

	      // Consider the incoming edge from the visited vertex
	      val minWeight = ghTree.edgesOf(nv).filter{ e =>
		val otherVertex = otherEndpoint(ghTree, nv, e)
		// Retain only those edges whose other vertex has been visited
		maxFlowMat.get(v.id, otherVertex.id) < Double.MaxValue
              }.map{ e =>
	        val ep = otherEndpoint(ghTree, nv, e)
	        min( maxFlowMat.get(v.id, ep.id), ghTree.getEdgeWeight(e) )
	      }.min

	      maxFlowMat.set(v.id, nv.id, minWeight)
	    }
	  }

	maxFlowMat.set( v.id, v.id, 0.0 )
	pb.update(n); n+=1
      }

    import java.io.{ FileWriter => JFileWriter }
    val ofile = new JFileWriter("MaxFlow.txt")

    for( i <- 0 until N ) {
      ofile.write( maxFlowMat.extractVector(true, i).getMatrix.getData.mkString(" ") + "\n" )
    }

    ofile.close

    maxFlowMat
  }

  def graphToDenseMatrix[V, E , G[V,E] <: WeightedPseudograph[V, E] ]( G : G[V,E]) = {
    val N = G.vertexSet.size
    val wmat = new SimpleMatrix(N,N)
    val vertIDMap = G.vertexSet.zipWithIndex.map{ vi => ( (vi._1, vi._2) )}.toMap

    // First, we build the weight matrix W and the diagonal D
    // Iterate over each vertex (row)
    G.vertexSet.foreach{
      vi =>
	var dweight = 0.0
        val cidx = vertIDMap(vi)
	// Iterate over v_i's adjacent edges
	G.edgesOf(vi).foreach{
	  eij =>
	    val (src, tgt, weight) = (G.getEdgeSource(eij), G.getEdgeTarget(eij), G.getEdgeWeight(eij))
	    val (i,j) = (vertIDMap(src), vertIDMap(tgt))
	    val pweight = 1.0 / G.degreeOf(vi)

	    if (i == cidx) {
	      wmat.set(i, j, pweight)
	    } else {
	      wmat.set(j, i, pweight)
	    }
	}

    }

    wmat
  }

  def graphToLaplacian[V, E , G[V,E] <: WeightedPseudograph[V, E] ]( G : G[V,E]) = {
    val N = G.vertexSet.size
    val wmat = new SimpleMatrix(N,N)
    val dmat = new SimpleMatrix(N,N)
    val vertIDMap = G.vertexSet.zipWithIndex.map{ vi => ( (vi._1, vi._2) )}.toMap

    // First, we build the weight matrix W and the diagonal D
    // Iterate over each vertex (row)
    G.vertexSet.foreach{
      vi =>
	var dweight = 0.0
        val cidx = vertIDMap(vi)
	// Iterate over v_i's adjacent edges
	G.edgesOf(vi).foreach{
	  eij =>
	    val (src, tgt, weight) = (G.getEdgeSource(eij), G.getEdgeTarget(eij), G.getEdgeWeight(eij))
	    val (i,j) = (vertIDMap(src), vertIDMap(tgt))
	    val pweight = 1.0 / G.degreeOf(vi)
	    dweight += pweight

	    if (i == cidx) {
	      wmat.set(i, j, pweight)
	    } else {
	      wmat.set(j, i, pweight)
	    }
	}

      dmat.set(cidx, cidx, dweight)
    }

    dmat.minus(wmat)
  }

  /*
     * Computes the all-pairs mean first passage time for the specified graph,
     * given an existing stationary probability distribution.
     * <P>
     * The mean first passage time from vertex v to vertex w is defined, for a
     * Markov network (in which the vertices represent states and the edge
     * weights represent state->state transition probabilities), as the expected
     * number of steps required to travel from v to w if the steps occur
     * according to the transition probabilities.
     * <P>
     * The stationary distribution is the fraction of time, in the limit as the
     * number of state transitions approaches infinity, that a given state will
     * have been visited. Equivalently, it is the probability that a given state
     * will be the current state after an arbitrarily large number of state
     * transitions.
     *
     * @param G
     *            the graph on which the MFPT will be calculated
     * @param edgeWeights
     *            the edge weights
     * @param stationaryDistribution
     *            the asymptotic state probabilities
     * @return the mean first passage time matrix
     */
    def computeMeanFirstPassageMatrix(
      G : WeightedPseudograph[SimpleVertex,SequenceWeightedEdge],
      weighted : Boolean ) =
	{

	  val N = G.vertexSet.size
	  val invNumEdges = 1.0 / G.edgeSet.size
	  val stationaryDistribution = G.vertexSet.toList.map{ v => G.degreeOf(v) * invNumEdges  }.toArray
	  println(N)
	  println(stationaryDistribution.size)
	  val temp = graphToDenseMatrix(G)

	  var i = 0
	  while ( i < temp.numRows ) {

	    var j = 0
	    while ( j < temp.numCols ) {
	      var value = -1 * temp.get(i,j) + stationaryDistribution(j)
	      if (i == j) { value += 1 }
	      if (value != 0) { temp.set(i,j,value) }
	      j +=1
	    }
	    i += 1

	  }

          val fundamentalMatrix = temp.invert
          temp.zero
	  i = 0;
	  while (i < temp.numRows) {

	     var j = 0
            while ( j < temp.numCols ) {
	      var value = -1.0 * fundamentalMatrix.get(i, j);
              value += fundamentalMatrix.get(j, j);
              if (i == j) { value += 1 }
              if (value != 0) { temp.set(i, j, value) }
	      j += 1
            }
	    i += 1
        }

        val stationaryMatrixDiagonal = new SimpleMatrix(N, N)
	(0 until N).foreach{  i => stationaryMatrixDiagonal.set(i, i, 1.0 / stationaryDistribution(i) ) }
	temp.mult(stationaryMatrixDiagonal)
    }



    def computeMyMeanFirstPassageMatrix(
      G : WeightedPseudograph[SimpleVertex,SequenceWeightedEdge],
      weighted : Boolean ) =
	{

	  val N = G.vertexSet.size
	  val vol = G.vertexSet.toList.map{ x => G.degreeOf(x).toDouble }.sum
	  val L = graphToLaplacian(G)
	  val weightedIdent = new SimpleMatrix( Array.fill(N,N){ 1.0 / N } )
          val fundamentalMatrix = L.minus(weightedIdent).invert.plus(weightedIdent)

	  weightedIdent.zero
	  var i = 0;
	  while (i < weightedIdent.numRows) {

	     var j = 0
            while ( j < weightedIdent.numCols ) {
	      weightedIdent.set(i,j,
				(fundamentalMatrix.get(i,i) + fundamentalMatrix.get(j,j) - 2 * fundamentalMatrix.get(i,j) ) )
	      j += 1
            }
	    i += 1
        }
	weightedIdent
    }

}
