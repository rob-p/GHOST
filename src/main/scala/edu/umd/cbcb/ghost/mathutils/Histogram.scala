package edu.umd.cbcb.ghost.mathtools

import scala.math._

object Histogram {
  implicit def toOption[T]( x:T ): Option[T] = { Option(x) }
}

class Histogram( vals: Array[Double], drange: Option[(Double,Double)] = None, nbins: Int = 10, normed: Boolean = false ) {

 // val testHist = new Histogram2(vals, drange, nbins)

  // The bin width.  Simply the bounds of the histogram divided by the number of bins
  private val h = {
    if ( !drange.isEmpty ) {
      (drange.get._2 - drange.get._1) / nbins
    } else {
      (vals.max - vals.min) / nbins
    }
  }

  val bins = (1 to nbins).map{ i => i * h }.toArray
  private val numVal = vals.size
  private val minVal = 0.0 //vals.min
  private val invNormFact = if (normed) { 1.0 / (h*numVal) } else { 1.0 }
  private val b = 0.0
  private def bin( v:Double ) = { min(max(0, floor((v - minVal - b) / h).toInt), nbins-1) }

  // Array of values for the histogram bins
  val hist = {
    val bv = Array.ofDim[Double](nbins)//DenseVector.zeros[Double](nbins)
    // For each value, compute the bin in which it resides and
    // appropriately update that bin's frequency (weight)
    vals.foreach{
      v =>
        if (v >=0 ){
          bv( bin(v) ) += invNormFact
        }
    }
    bv
  }
}
