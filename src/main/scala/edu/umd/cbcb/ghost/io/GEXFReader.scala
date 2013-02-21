package edu.umd.cbcb.ghost.io

import scala.collection.mutable._
import scala.math._
import scala.xml._

import org.jgrapht.graph._
import edu.umd.cbcb.ghost.graphtools.SequenceWeightedEdge

import java.io.FileNotFoundException

class GEXFReader(filename: String, useWeight: Boolean = true) {
  type GraphT = WeightedPseudograph[SimpleVertex, SequenceWeightedEdge]
  val G = new GraphT( classOf[ SequenceWeightedEdge ] )
  val nameVertexMap = new HashMap[String, SimpleVertex]
  var attributeIDMap = scala.collection.immutable.Map.empty[String, Int]

  def parse: GraphT = {

    var tlSeq: Elem = null

    try { 
      tlSeq = XML.loadFile(filename)
    } catch { 
      case nsfe: FileNotFoundException => println("No Such File Exception; reason = %s".format(nsfe.getMessage) )
    }

    val graphSeq = tlSeq \ "graph"
    val attrSeq = graphSeq \\ "attribute"

    attributeIDMap = attrSeq.map{ x => (x.attribute("title").get.text, x.attribute("id").get.text.toInt) }.toMap[String, Int]
    
    val nodeSeq = graphSeq \\ "node"
    
    nodeSeq.zipWithIndex.foreach{ 
      xi => 
	val label = xi._1.attribute("label").get.text
	val id = xi._2
        val attSeq = xi._1 \\ "attvalue"
	val gnameAttr = attSeq.filter{ a => a.attribute("for").get.text.toInt == attributeIDMap("gname") }
        val gname = gnameAttr(0).attribute("value").get.text
        val vert = new SimpleVertex(gname, id)
        nameVertexMap.getOrElseUpdate(label, vert)
        G.addVertex(vert)
    }
    
    val edgeSeq = graphSeq \\ "edge"
    val defaultNode = new Text("1.0")
    edgeSeq.foreach{ 
      e =>
	val src = nameVertexMap( e.attribute("source").get.text )
        val tgt = nameVertexMap( e.attribute("target").get.text )
	val ne = new SequenceWeightedEdge
        val added = G.addEdge( src, tgt, ne )
        val watt = e.attribute("weight")
        val weight = if ( watt.isEmpty || (! useWeight) ) { 1.0 } else {  watt.get.text.toDouble }
        G.setEdgeWeight( ne, weight )
    }

    G
  }

}

