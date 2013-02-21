package edu.umd.cbcb.ghost.graphtools

import org.jgrapht.graph.DefaultWeightedEdge

class SequenceWeightedEdge( _seqWeight: Double = 0.0 ) extends DefaultWeightedEdge { 
  var seqWeight = _seqWeight
}
