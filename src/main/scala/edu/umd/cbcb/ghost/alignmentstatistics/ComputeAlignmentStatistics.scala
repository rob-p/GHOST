package edu.umd.cbcb.ghost.alignmentstatistics

import scopt._
import scala.io._
import scalax.file.Path
import scalax.file.ImplicitConversions._
import scala.math._
import scala.xml.XML._
import scala.collection.JavaConversions._
import scala.collection.mutable.{ Map => MMap }
import org.jgrapht.graph._
import java.{ io => jio }
import java.io.{ FileWriter => JFileWriter, BufferedWriter => JBufferedWriter, PrintWriter => JPrintWriter }

import edu.umd.cbcb.ghost.io.{ GEXFReader, EdgeListReader, SimpleVertex }
import edu.umd.cbcb.ghost.graphtools.{ Utils, SequenceWeightedEdge }
import edu.umd.cbcb.ghost.align.{ Alignment }

object ComputeAlignmentStatistics { 
  type VertexT = SimpleVertex
  type EdgeT = SequenceWeightedEdge
  type StringEdgeT = Edge
  type GraphT = WeightedPseudograph[SimpleVertex, SequenceWeightedEdge]
  type SubgraphT = Subgraph[SimpleVertex, SequenceWeightedEdge, GraphT] 

  /**
   * Configuration Options
   */
  class Config { 
    var alignment = ""
    var g1FileName = ""
    var g2FileName = ""
    var alignerName = ""
  }
  
  /**
   * Give ourselves a nice unordered edge
   */
  class Edge( _u: String, _v: String ) { 
    val (u,v) = if ( _u < _v ) { (_u,_v) } else { (_v, _u) }

    // We should be able to hash edges
    override def hashCode = { (u,v).hashCode }
    // And compare them for equality
    override def equals(o: Any) = {
      if (o.isInstanceOf[Edge]) { 
	val oe = o.asInstanceOf[Edge]
	oe.u == u && oe.v == v
      } else { false }
    }
    // And print them out nicely
    override def toString(): String = { 
      "Edge[%s]".format((u,v))
    }
  }

  def main( args: Array[String] ) = { 
    var config = new Config
    val parser = new OptionParser("ComputeAlignmentStatistics") { 
      opt("u", "g1", "first input graph",
	{ v: String => config.g1FileName = v})
      opt("v", "g2", "second input graph",
	{ v: String => config.g2FileName = v})
      opt("a", "alignment", "alignment file or directory",
	{ v: String => config.alignment = v})
      opt("n", "name", "aligner name",
        { v: String => config.alignerName = v})
    }

    if (parser.parse(args)) {
      
      val alignmentPath = string2path(config.alignment)
      var alignmentFiles = List.empty[Path]

      if (alignmentPath.isDirectory) { 
        alignmentFiles = alignmentPath.descendants().filter{ x => x.name.endsWith(".aln") } ++: alignmentFiles
      } else {
        alignmentFiles = alignmentPath :: alignmentFiles
      }


      System.err.println("<alignments method=\""+config.alignerName+"\">")
      for ( afile <- alignmentFiles ) { 
        // Read the alignment file
        val ws = """\s+""".r
        val alnStr = Source.fromFile( afile.path ).getLines.map{ x => val Array(l,r) = ws.split(x); (l,r) }.toMap

        // Create a new GEXFReader and read in the graphs
        val g1r = new GEXFReader(config.g1FileName, false)
        val g1f = g1r.parse
        val g2r = new GEXFReader(config.g2FileName, false)
        val g2f = g2r.parse
        
        // Remove self-loops from the graphs
        /*
        val remE1 = g1f.edgeSet.filter{ e => g1f.getEdgeSource(e) == g1f.getEdgeTarget(e) }
        val remE2 = g2f.edgeSet.filter{ e => g2f.getEdgeSource(e) == g2f.getEdgeTarget(e) }
        g1f.removeAllEdges( remE1 )
        g2f.removeAllEdges( remE2 )
        */

        var alignment = Alignment.empty[SimpleVertex]

        // f : G => H
        val fwMap = alnStr.map{ case (x,y) => 
          val (u,v) = (g1r.nameVertexMap.get(x), g2r.nameVertexMap.get(y))
          (u,v) }.filterNot{ case (x,y) => x == None || y == None }.map{ x => 
            alignment.insert(x._1.get,x._2.get,0.0)
            (x._1.get, x._2.get) 
          }.toMap
        // f^{-1} : H => G
        val revMap = fwMap.map{ x => (x._2, x._1) }.toMap
      
      val vset1: Set[VertexT] = fwMap.keys.toSet
      val vset2: Set[VertexT] = fwMap.values.toSet
      var sg1 = new SubgraphT(g1f, vset1 ) // Should just be G
      var sg2 = new SubgraphT(g2f, vset2 ) // Ind( f(G) ) 

      var mappedEdges = sg1.edgeSet.map{ e => new Edge(fwMap(g1f.getEdgeSource(e)).name, fwMap(g1f.getEdgeTarget(e)).name) }.toSet
      var inducedEdges = sg2.edgeSet.map{ e => new Edge(g2f.getEdgeSource(e).name, g2f.getEdgeTarget(e).name) }.toSet
      val sharedEdges = mappedEdges & inducedEdges
      // val unsharedEdges = inducedEdges &~ mappedEdges

      val inducedSubgraph = new GraphT( classOf[EdgeT] )
      val nameVertMap = MMap.empty[String, SimpleVertex]
      sg2.vertexSet.foreach{ v => inducedSubgraph.addVertex(v); nameVertMap(v.name) = v }
      inducedEdges.foreach{ e => val ne = new EdgeT(); inducedSubgraph.addEdge( nameVertMap(e.u), nameVertMap(e.v), ne ) } 
      
      val commonSubgraph = new GraphT( classOf[EdgeT] )
      sg2.vertexSet.foreach{ v => commonSubgraph.addVertex(v) }
      sharedEdges.foreach{ e => val ne = new EdgeT(); commonSubgraph.addEdge( nameVertMap(e.u), nameVertMap(e.v), ne ) } 
      val lcc = Utils.largestConnectedComponent(commonSubgraph)
      val lccInduced = new SubgraphT(g2f, lcc.vertexSet)

      val gmap = new DefaultGraphMapping(fwMap, revMap, sg1, sg2)

      val fmapped = sg1.edgeSet.count{ e1 => gmap.getEdgeCorrespondence(e1,true) != null }
      val (s1, s2) = (sg1.edgeSet.size, sg2.edgeSet.size)

      val ec = fmapped / s1.toFloat
      val aln = 
      <alignment>
        <alpha>{ afile.path }</alpha>
        <ec>{ fmapped }</ec>
        <ics>{ Utils.inducedConservedStructure(alignment.->-!,g1f,g2f) * 100.0 }</ics>
        <lccs>
          <vertices>{ lcc.vertexSet.size }</vertices>
          <edges>{ lcc.edgeSet.size }</edges>
        </lccs>
      </alignment>;
      System.err.println(aln.toString())

      /*
      println("===== Alignment "+afile.path+" ======")
      println("One Way Edge Accuracy is %f %% (%d of %d edges)".format( ec * 100.0, fmapped, s1 ))
      println("Symmetric Edge Accuracy is %f %%".format( ec * (fmapped.toFloat / inducedSubgraph.edgeSet.size) * 100.0))
      println("Symmetric Edge Accuracy of LCCS %f %%".format( (lcc.edgeSet.size.toFloat / lccInduced.edgeSet.size) * 100.0 ) )
      println("LCCS has %d vertices and %d edges edges".format(lcc.vertexSet.size, lcc.edgeSet.size) )
      println("ICS Score %f".format( (fmapped.toFloat / inducedSubgraph.edgeSet.size)))
      */
      //println("EC is %f".format( ec ))
      }
      System.err.println("</alignments>")
    }
  }

}
