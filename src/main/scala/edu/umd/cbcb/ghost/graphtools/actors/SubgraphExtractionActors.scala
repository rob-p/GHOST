package edu.umd.cbcb.ghost.graphtools.actors

import scala.collection.mutable.{ HashSet => MHashSet, OpenHashMap => MHashMap, MultiMap }
import scala.collection.immutable.{ TreeMap }
import scala.collection.JavaConversions._
import scala.collection.parallel.mutable.ParArray

// import scala.actors.Actor
// import scala.actors._

import akka.actor.{ Actor, ActorRef }
import akka.actor.{ DeadLetter }
// Akka 1.x
//import akka.dispatch.{ MessageDispatcher }
//import akka.event.EventHandler

import scala.math._

import org.jgrapht.graph._
import org.jgrapht.traverse._
import org.ejml.simple._

import java.io.{ File => JFile, FileWriter => JFileWriter }
import java.io.{ FileOutputStream => JFileOutputStream, FileInputStream => JFileInputStream }
import java.io.{ BufferedOutputStream => JBufferedOutputStream }
import java.io.{ DataOutputStream => JDataOutputStream, DataInputStream => JDataInputStream, PrintStream => JPrintStream }

import net.robpatro.utils.console.ProgressBar
import edu.umd.cbcb.ghost.io.{ SimpleVertex }
import edu.umd.cbcb.ghost.graphtools.{ LaplacianCalculator, SequenceWeightedEdge }

/* Messages accepted by SubgraphExtractor and SubgraphExtractorReceiver */

// Tells an extractor to extract the subgraph centered at v and return the descriptor
case class ExtractSubgraphDescriptor(v: SimpleVertex)

// Companion object of SubgraphDescriptor
object SubgraphDescriptor {
  def readFromInputStream(descFile: JDataInputStream, maxGlobalLevel: Int, nodeId: Int = 0) = {
    // write the vertex name
    val nvert = new SimpleVertex(descFile.readUTF(), nodeId)
    val maxLevel = maxGlobalLevel

    val vertexSets = new MHashMap[Int, collection.Set[SimpleVertex]]
    val spectra = new MHashMap[Int, Array[Double]]
    val densities = new MHashMap[Int, Double]

    for (level <- 1 to maxLevel) {
      // read the level
      val level = descFile.readInt

      // read the vertex set
      val numVerts = descFile.readInt()
      val vset = new MHashSet[SimpleVertex]
      (0 until numVerts).foreach { i => vset.add(new SimpleVertex(descFile.readUTF(), i)) }
      vertexSets(level) = vset

      // read the spectrum
      val specSize = descFile.readInt
      //val spectrum = Array.ofDim[Double](specSize)
      //spectra(level) = spectrum
      spectra(level) = (0 until specSize).map { i => descFile.readDouble }.toArray

      // write the density
      densities(level) = descFile.readDouble
    }
    //print("Vertex "+nvert.name+" spectral sums = ("+spectra.map{ kv => kv._2.sum }.mkString(",")+")")
    new SubgraphDescriptor(nvert, vertexSets, spectra, densities, maxLevel)
  }

  def readFromInputStream(descFile: JDataInputStream): (ParArray[SubgraphDescriptor], Int) = {
    val numDesc = descFile.readInt
    val maxLevel = descFile.readInt
    val descs = (0 until numDesc).map { i => SubgraphDescriptor.readFromInputStream(descFile, maxLevel, nodeId = i) }.toArray.par
    (descs, maxLevel)
  }

}

// The descriptor returned from the message aboves
class SubgraphDescriptor(val centerVertex: SimpleVertex,
  val vertexSets: collection.mutable.Map[Int, collection.Set[SimpleVertex]] = new MHashMap[Int, collection.Set[SimpleVertex]],
  var spectra: MHashMap[Int, Array[Double]] = new MHashMap[Int, Array[Double]],
  val densities: MHashMap[Int, Double] = new MHashMap[Int, Double],
  var maxLevel: Int = 0) {

  def addLevel(l: Int, vs: scala.collection.Set[SimpleVertex], evals: Array[Double], density: Double) = {
    vertexSets(l) = vs
    spectra(l) = evals
    densities(l) = density
    maxLevel = max(l, maxLevel)
  }

  def copyPreviousLevel(l: Int) = {
    val prevLevel = l - 1
    vertexSets(l) = vertexSets(prevLevel)
    spectra(l) = spectra(prevLevel)
    densities(l) = densities(prevLevel)
    maxLevel = max(l, maxLevel)
  }

  def writeToSimpleFile(descFile: JPrintStream) = {
    // write the vertex name
    descFile.println(centerVertex.name)
    var totVals = 0
    for (level <- 1 to maxLevel) {
      totVals += vertexSets(level).size
      descFile.println(totVals.toString)
      descFile.println(spectra(level).mkString(" "))
    }

  }

  def writeToOutputStream(descFile: JDataOutputStream, dnum: Int) = {

    // write the vertex name
    descFile.writeUTF(centerVertex.name)

    for (level <- 1 to maxLevel) {
      // write the level
      descFile.writeInt(level)

      // write the vertex set
      descFile.writeInt(vertexSets(level).size)
      vertexSets(level).foreach { v => descFile.writeUTF(v.name) }

      // write the spectrum
      descFile.writeInt(spectra(level).size)
      spectra(level).foreach { v => descFile.writeDouble(v) }

      // write the density
      descFile.writeDouble(densities(level))
    }
    descFile.flush()

  }

}

case class SubgraphDescriptorMsg(d: SubgraphDescriptor)

// Query if the actor is done
case class CompleteQuery()
// Tells the actor that they are done
case class Stop()

class SubgraphDescriptorReceiver(reqNumDesc: Int, k: Int, descFileStream: Either[JDataOutputStream, JPrintStream]) extends Actor {

  import context._

  val simpleOutput = descFileStream.isRight
  var numReceivedDesc = 0
  val pb = new ProgressBar(reqNumDesc, "*")

  // Write the total number of descriptors to the top of the file
  if (simpleOutput) {
    descFileStream.right.get.println(reqNumDesc.toString)
    descFileStream.right.get.println(k.toString)
  } else {
    descFileStream.left.get.writeInt(reqNumDesc)
    descFileStream.left.get.writeInt(k)
  }

  def receive = {
    case SubgraphDescriptorMsg(desc) => {
      writeDescriptor(desc)
      pb.update(numReceivedDesc)
      if (numReceivedDesc == reqNumDesc) { context.stop(self) }
    }
    case CompleteQuery() => {
      // Akka 1.x
      //self.reply(done)
      sender ! done
    }
    case Stop() => {
      pb.done()
    }
    case d: DeadLetter => println(d)
    case _ => println("received unkown message")
  }

  /*
  def act() = {

    loopWhile( !done ) {
      react {
	case SubgraphDescriptorMsg( desc ) => {
	  writeDescriptor(desc)
	  pb.update(numReceivedDesc)
	}
	case Stop() => {
	  pb.done()
	  exit()
	}
      } // end react
    } // end loopWhile

  }
  */
  def done() = { numReceivedDesc == reqNumDesc }

  def writeDescriptor(desc: SubgraphDescriptor) = {
    if (simpleOutput) {
      desc.writeToSimpleFile(descFileStream.right.get)
    } else {
      desc.writeToOutputStream(descFileStream.left.get, numReceivedDesc)
    }
    numReceivedDesc += 1
  }

}

class SubgraphExtractor[T <: AbstractGraph[SimpleVertex, SequenceWeightedEdge]](
  _G: T,
  _k: Double,
  _rec: ActorRef,
  _mfpt: Option[SimpleMatrix]) extends Actor {
  // Akka 1.x
  //_dispatcher : MessageDispatcher)  extends Actor {
  /*
   * The graph and k hop neighborhood with which
   * this actor will work
   */
  val G = _G
  val k = _k
  val rec = _rec
  val mfpt = _mfpt
  // Akka 1.x
  //self.dispatcher = _dispatcher

  def receive = {
    case ExtractSubgraphDescriptor(vert) => rec ! SubgraphDescriptorMsg(generateSubgraphDescriptor(vert))
    case d: DeadLetter => println(d)
    case _ => println("received unknown message")
  }
  /*
   def act() = {
    loop {
      react {
	case ExtractSubgraphDescriptor( vert ) => rec ! SubgraphDescriptorMsg( generateSubgraphDescriptor( vert ) )
	case Stop => exit()
      }
    }
  }
  */

  def generateSubgraphDescriptor(vert: SimpleVertex): SubgraphDescriptor = {
    type SubgraphT = Subgraph[SimpleVertex, SequenceWeightedEdge, T]
    val desc = new SubgraphDescriptor(vert)
    val cfi = new ClosestFirstIterator(G, vert, k)
    var vertexSet = new MHashSet[SimpleVertex]
    var distVertMap: Map[Int, Set[SimpleVertex]] = cfi.toSet.groupBy { i => cfi.getShortestPathLength(i).toInt }
    // Ensure that the distVertMap has entries, even if empty, for all levels
    for (chop <- 1 to k.toInt; if !(distVertMap contains chop)) {
      distVertMap += (chop -> Set.empty[SimpleVertex])
    }

    val useMfpt = !mfpt.isEmpty

    vertexSet.add(vert)
    for (chop <- 1 to k.toInt) {
      var newVerts, levelVerts = distVertMap(chop).toList
      // If we're using the max-flow criterion for determining neighbors
      if (useMfpt) {
        val mfptVec = mfpt.get.extractVector(true, vert.id)
        levelVerts = levelVerts.filterNot { x => mfptVec.get(x.id).isNaN }
        val sfact = 1.0 / levelVerts.size
        val mean = levelVerts.map { x => mfptVec.get(x.id) }.sum * sfact
        val maxFlow = levelVerts.map { x => mfptVec.get(x.id) }.max.toDouble
        val stddev = sqrt(levelVerts.map { x => sfact * pow(mfptVec.get(x.id) - mean, 2.0) }.sum)
        val cutoff = mean - 3.0 * stddev
        newVerts = levelVerts.filter { x => mfptVec.get(x.id) >= cutoff }
      }

      val psize = vertexSet.size
      vertexSet ++= newVerts
      if ( chop == 1 || vertexSet.size > psize ) { 
        val N = vertexSet.size
        val subGraph = new SubgraphT(G, vertexSet)
        val density = N.toFloat / subGraph.edgeSet.size

        val lapCalc = new LaplacianCalculator(subGraph, (x, y) => x)
        val (lap, vols) = lapCalc.compute

        /* Use the dense eigensolver */
        import org.ejml.ops._
        import org.ejml.data.DenseMatrix64F
        import org.ejml.alg.dense.decomposition._

        val ofile = new java.io.BufferedWriter(new java.io.FileWriter("EigenDecompTest/"+vert.name+"_"+chop+".txt"))
        var crow = 0
        
        val denseL = new DenseMatrix64F(N, N)
        val it = lap.iterator
        while (it.hasNext) {
          val ent = it.next
          denseL.unsafe_set(ent.row, ent.column, ent.get)
          if (crow != ent.row) { ofile.write("\n"); crow = ent.row }
          ofile.write(ent.column+":"+ent.get+" ") 
        }

        ofile.close();
        /*
      val rowIt = lap.iterator
      while( rowIt.hasNext ){
	val row = rowIt.next
	val i = row.index
	val colIt = row.vector.iterateNonZero
	while( colIt.hasNext ){
	  val elem = colIt.next
	  val j = elem.index
	  denseL.unsafe_set(i, j, elem.get)
	}
      }
      */
        // Perform the eigendecomposition
        val decomposer = DecompositionFactory.eigSymm(N, true)
        decomposer.decompose(denseL)

        // Sort the indices
        val sortedInds = (0 until N).map { i => (decomposer.getEigenvalue(i).real, i) }.toArray.sorted.map { vi => vi._2 }
        // Get the eigenvalues and eigenvectors in the sorted order
        val devals = (0 until N).map { i => decomposer.getEigenvalue(sortedInds(i)).real }.toArray
        val sigvals = devals
        // val sigvals  = if (N > 1) { decomposer.getEigenVector( sortedInds(1) ).getData.toArray } else { Array.ofDim[Double](1) }
        /*(1 until chop).map{
	  i =>
	    val a : Array[Double] = if( i < N){
	      decomposer.getEigenVector( sortedInds(i) ).getData.toArray
	    } else {
	      Array.ofDim[Double](N)
	    }
            (i,a)
        }.toMap
      */

        // Add this information to the vertex's discriptor
        desc.addLevel(chop, distVertMap(chop), sigvals, density)
      } else {
        desc.copyPreviousLevel(chop)
      }
      /*sigvals.keys.foreach{
	i =>
	  desc.addLevel(i, vertexSet, sigvals(i), density)
      }*/

    }

    desc
  }

}
