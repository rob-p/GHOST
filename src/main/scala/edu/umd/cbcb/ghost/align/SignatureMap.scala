package edu.umd.cbcb.ghost.align

import edu.umd.cbcb.ghost.graphtools.actors._
import edu.umd.cbcb.ghost.mathtools.Histogram
import edu.umd.cbcb.ghost.mathutils.UtilFunctions

import scala.collection.JavaConversions._
import scala.collection.parallel.mutable.ParArray
import scala.collection.mutable.{ ArrayBuffer, HashMap => MHashMap, HashSet => MHashSet }
import scala.math._
import java.util.concurrent.{ ConcurrentHashMap => CHMap }

import java.io.{ DataInputStream => JDataInputStream, BufferedInputStream => JBufferedInputStream, FileInputStream => JFileInputStream}
import java.util.zip.{ GZIPInputStream => JGZIPInputStream }

import org.jgrapht.graph._
import org.jgrapht.traverse._
import edu.umd.cbcb.ghost.io.{ SimpleVertex }
import edu.umd.cbcb.ghost.graphtools.{ LaplacianCalculator, SequenceWeightedEdge }

import org.ejml.data.DenseMatrix64F
import org.ejml.simple.SimpleMatrix

class SignatureMap( val name: String, val G: WeightedPseudograph[SimpleVertex, SequenceWeightedEdge], val numBins: Int ) {
  type GraphT = WeightedPseudograph[SimpleVertex, SequenceWeightedEdge]
  type SubgraphT = Subgraph[SimpleVertex, SequenceWeightedEdge, GraphT]

  val nameVertexMap = G.vertexSet.map{ v => (v.name, v) }.toMap
  var dmat = new SimpleMatrix(0,0)
  var nameToId: collection.concurrent.Map[String, Int] = new CHMap[String, Int]()
  var idToName: collection.concurrent.Map[Int, String] = new CHMap[Int, String]()
  var avgDensity: collection.concurrent.Map[String, Double] = new CHMap[String, Double]()
  var spectra: collection.concurrent.Map[String, ArrayBuffer[Array[Double]]] = new CHMap[String, ArrayBuffer[Array[Double]]]()
  var size = 0
  var binWidths = Array.empty[Int]
  var maxHop = 0
  var maxDeg = G.degreeOf( G.vertexSet.maxBy{ x => G.degreeOf(x) } ) + 1
  val avgDeg = G.vertexSet.toList.map{ x => G.degreeOf(x) }.sum.toFloat / G.vertexSet.size

  private var neighborhoods = new CHMap[SimpleVertex, CHMap[Int, scala.collection.Set[SimpleVertex]] ]()

  def readSignatures( sigFileName: String, histBase: Double ) = {
    // Read the signatures from file
    val sigFile = new JDataInputStream(new JGZIPInputStream( new JBufferedInputStream( new JFileInputStream( sigFileName ) ) ) )
    var (sigs, khop) = SubgraphDescriptor.readFromInputStream(sigFile)
    sigFile.close()
    maxHop = khop
    // Compute the histograms
    var (miSigSum, maSigSum) = (Double.MaxValue,-Double.MaxValue)

    binWidths = (1 to khop).map{ h =>  h*max(numBins, 2*avgDeg.toInt) }.toArray //max(40, (pow(histBase,h) * 0.1).toInt) }.toArray// }max(100, pow(histBase,h).toInt) }.toArray
    val ncols = binWidths.sum
    val nrows = sigs.size
    val dat = Array.ofDim[Double](nrows, ncols)
    sigs.zipWithIndex.foreach{
      si =>
      val (s,i) = si
      val (name, id) = (s.centerVertex.name, s.centerVertex.id);
      nameToId += (name -> id)
      idToName += (id -> name)

      spectra(name) = ArrayBuffer.empty[Array[Double]]
      val gamma = 0.005

      neighborhoods( nameVertexMap(name) ) = new CHMap[Int, collection.Set[SimpleVertex]]()
      val specBuf = ArrayBuffer.empty[Double]
      val lvlNorm = 1.0 / s.maxLevel
      val ra = (1 to s.maxLevel).foreach{
        level =>
          neighborhoods( nameVertexMap(name) )( level ) = s.vertexSets(level).map{ v => nameVertexMap(v.name) }
          val spectrum = s.spectra(level)
          specBuf ++= spectrum
          val lgamma =  gamma // level
          spectra(name) += UtilFunctions.ipsenMikhailovVector(spectrum, gamma=lgamma)
          val density = s.densities(level)
          if (level == 1) { avgDensity += (s.centerVertex.name -> lvlNorm * density) } else { avgDensity(s.centerVertex.name) += lvlNorm * density }
          //val hist = new Histogram(spectrum, drange=Option((0.0,2.0)), nbins=binWidths(level-1) )
          //val scaleFact = 1.0 / sqrt( hist.hist.map{ x => x*x }.sum )
          //val r = hist.hist.map{ h => h.toDouble * scaleFact } // * (0.25 / density) }
          //hist.hist
          //r
          //}.toArray
      }
      //spectra(name) = UtilFunctions.ipsenMikhailovVector(specBuf.toArray)
      //dat(i) = ra
    }

    //dmat = new SimpleMatrix( dat )
    size = nrows
  }

  def getName(id: Int) = { idToName(id) }
  def getId(name: String) = { nameToId(name) }

  def neighborhoodByName(name: String, level: Int) = {
    val cvertex = nameVertexMap(name)
    val cvertexMap = neighborhoods(cvertex)
    var vset = Set( cvertex )
    (1 to level).foreach{ l => vset |= cvertexMap(l) }
    vset
  }

  def removeFromNeighborhoods( rset: Set[Int] ) = {
    var toRemove = new MHashSet[(SimpleVertex, Int)]
    rset.foreach{
      vi =>

      val cv = nameVertexMap( idToName(vi) )
      neighborhoods.foreach{
        kv =>
          kv._2.foreach{
                lset =>
                  if (lset._2 contains cv) {
                      toRemove += ((kv._1, lset._1))
                  }
          }
      }
      toRemove.foreach{ kl => neighborhoods(kl._1)(kl._2) -= cv }
      if (neighborhoods contains cv) { neighborhoods -= cv }

    }
  }

}

