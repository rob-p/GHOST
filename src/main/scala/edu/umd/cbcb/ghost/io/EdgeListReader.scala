package edu.umd.cbcb.ghost.io

import scala.collection.mutable._
import scala.math._
import org.jgrapht.graph._

class EdgeListReader(filename: String) {

  val G = new SimpleWeightedGraph[SimpleVertex, DefaultWeightedEdge]( classOf[ DefaultWeightedEdge ] )
  val nameVertexMap = new HashMap[String, SimpleVertex]

  def readEdge(lineit: Iterator[String]): Option[(SimpleVertex, SimpleVertex, Double)]  = {

    if ( lineit.hasNext ) { 
      val line = lineit.next
    // The attrbute dictionary is a python-esque dictionary
    // { 'key0' : val, 'key1' : val }
    val hasAttriuteDict = line.contains('{')
    val attributeDictIndex = if ( hasAttriuteDict ) { line.indexOf('{') } else { line.length }
    val st = line.substring(0,attributeDictIndex).split(' ').map{ x => x }
    val (src, tgt) = (st(0), st(1))
    var weight = 1.0

    if ( hasAttriuteDict ) { 
      val attributeDictList = line.substring(attributeDictIndex+1, line.length-1).split(", ").toList
      val attribMap = attributeDictList.map{ 
	kv =>
	  val toks = kv.split(": ")
	  ( toks(0).substring(1,toks(0).length-1), toks(1) )
      }.toMap
      weight = attribMap("weight").toDouble
    }

    val srcVertex = nameVertexMap.getOrElseUpdate(src, new SimpleVertex(src, nameVertexMap.size))
    val tgtVertex = nameVertexMap.getOrElseUpdate(tgt, new SimpleVertex(tgt, nameVertexMap.size))

    Some(srcVertex, tgtVertex, weight)

    } else { 
      None
    }

  }

  def parse: SimpleWeightedGraph[SimpleVertex, DefaultWeightedEdge]  = {

    val fin = io.Source.fromFile(filename)
    var lineit = fin.getLines()
    var nentry : Option[(SimpleVertex, SimpleVertex, Double)] = None

    while ( {nentry = readEdge(lineit); nentry != None} ) {
      val (src, tgt, weight) = nentry.get
      if ( !G.containsVertex(src) ) { G.addVertex(src) }
      if ( !G.containsVertex(tgt) ) { G.addVertex(tgt) }
      val e = G.addEdge(src, tgt)
      G.setEdgeWeight(e, weight)
    }
    
    G
  }

}

