package edu.umd.cbcb.ghost.alignmentstatistics

import scala.io._
import scala.math._

object CompareExactAlignments{ 

  def main( args: Array[String] ) = { 
    val space = "\\s+".r
    val alignment0 = Source.fromFile(args(0)).getLines.map{ l => val Array(i,j) = space.split(l); (i,j) }.toSet
    val alignment1 = Source.fromFile(args(1)).getLines.map{ l => val Array(i,j) = space.split(l); (i,j) }.toSet
    val shared = alignment0 & alignment1
    val ms = min(alignment0.size, alignment1.size)
    println("Shared Assignments between %s and %s : %d = %f %%".format(args(0), args(1), shared.size, 100.0 * (shared.size/ms.toFloat)) )
    println( shared.mkString(", ") )
  }

}
