package edu.umd.cbcb.ghost.mathutils

import scala.collection.mutable.{ ArrayBuffer }
import scala.math._

import org.ejml.simple._
import org.ejml.data._
import org.ejml.ops.NormOps

object UtilFunctions {

  def cosineDistance( v1 : SimpleMatrix, v2 : SimpleMatrix ) = {
    1.0 - (v1.dot(v2) / ( NormOps.normP2(v1.getMatrix) * NormOps.normP2(v2.getMatrix) ))
  }

  def l1Distance( wvec : SimpleMatrix, v1: SimpleMatrix, v2 : SimpleMatrix ) = {
    (NormOps.normP1( v1.minus(v2).elementMult(wvec).getMatrix ))
  }

  def klDiv( p1:Array[Double], p2:Array[Double] ) = {
    var s = 0.0
    var i = 0
    val N = p1.size
    while (i < N) {
      if ( p1(i) > 0 && p2(i) != 0.0 ) { 
        s += p1(i) * log( p1(i) / p2(i) )
      }
      i += 1
    }
    max( 0.0, s )
  }

  def jsDist( p1:Array[Double], p2:Array[Double] ) = {

    def avgDist( d1:Array[Double], d2:Array[Double] ) = {
      var i = 0
      val N = d1.size
      val d = Array.fill(d1.size)(0.0)
      while (i < N) { d(i) = 0.5 * (d1(i)+d2(i)); i += 1 }
      d
    }

    val p = avgDist(p1, p2)
    val (k1, k2) = ( 0.5 * klDiv(p1,p), 0.5 * klDiv(p2,p) )
    sqrt(k1 + k2)
  }

  def structDistance( imVec1:ArrayBuffer[Array[Double]], imVec2:ArrayBuffer[Array[Double]] ) = {
    var dist = 0.0
    var kthdist = 0.0
    var j = 0
    val M = imVec1.size
    while( j < M ) {
      val (v1, v2) = (imVec1(j), imVec2(j))
      val d = jsDist(v1,v2)
      dist += d
      j += 1
    }
    assert(dist < M)
    dist / M
  }

  // Compute the density (i.e. \rho) of the spectrum w.r.t. the given omega
  def lorentzDensity( spec:Array[Double], omega:Double, gamma:Double ) = {
    var d = 0.0
    val gammasq = (gamma*gamma)
    spec.foreach{ omegai =>
      d += (gamma)/( ((omega-omegai)*(omega-omegai)) + gammasq )
    }
    d
   }

  def gaussianDensity( spec:Array[Double], omega:Double, gamma:Double ) = { 
    val sigma = gamma
    val sigmasq = sigma * sigma
    // Assume a minimum density to avoid divide by zero errors
    var d = 1e-30
    val norm = 1.0 / sqrt(2.0*math.Pi*sigmasq)
    spec.foreach{ omegai =>
      d += norm * math.exp(-((omega-omegai)*(omega-omegai)) / (2.0*sigmasq))
    }
    d
  }

  def ipsenMikhailovVector( spec: Array[Double], gamma:Double=0.2 ) = {
    val delta = 0.1*gamma
    val ub = min(201.0,(200 + delta * 100.0)).toInt
    val step = max(1.0, (delta * 100.0)).toInt
    val omegas = Range(0, ub, step).map{ _/100.0 }
    var dens = Array.fill(omegas.size)({0.0})
    var i = 0
    var c = 0.0
    val N = omegas.size
    while( i < N ){
      val omega = omegas(i)
      //val d = lorentzDensity(spec, omega, gamma)
      val d = gaussianDensity(spec, omega, gamma)
      c += d
      dens(i) = d
      i += 1
    }
    val k = 1.0/c
    dens = dens.map{ _*k }
    dens
  }

  def ipsenMikhailovDistance( imVec1:ArrayBuffer[Array[Double]], imVec2:ArrayBuffer[Array[Double]] ) = {
    // Compute the distance
    var dist = 0.0
    var kthdist = 0.0
    var j = 0
    val M = imVec1.size
    while( j < M ) {

      val (v1, v2) = (imVec1(j), imVec2(j))
      var i = 0
      val N = imVec1.size
      while( i < N ) {
        val d = v1(i) - v2(i)
        kthdist += d*d
        i += 1
      }
      j += 1
      dist += sqrt(kthdist)

    }
    dist / M
  }

}



/*
def quadraticChi( l: SimpleMatrix, r: SimpleMatrix, A: SimpleMatrix, m: Double ) = {
  val Z = l.plus(r).mult(A)
  val N = Z.numRows
  val M = Z.numCols

  // Z ( Z == 0 ) = 1 && Z.^m && Z = 1.0 ./ Z
  for( i <- 0 until N; j <- 0 until M ) {
    val v = if (Z.get(i,j) == 0.0) { 1.0 } else { Z.get(i,j) }
    Z.set(i, j, 1.0 / pow(v,m) )
  }

  // (l-r) * (1.0 / Z)
  val D = (l.minus(r)).elementMult(Z)

  // dist = sqrt( max( D * A * D', 0.0 ) )
  sqrt( max( D.mult(A).mult(D.transpose()).get(0,0), 0.0) )
}

def  klDiv( l: Array[Double], r: Array[Double]) = {
  (0 until l.size).foldLeft(0.0){ (total, i) =>
    if (r(i) > 0 && l(i) > 0) { total + (l(i) * log( l(i) / r(i) )) } else { total }
			       }
}
*/
