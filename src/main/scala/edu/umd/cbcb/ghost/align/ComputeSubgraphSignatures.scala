package edu.umd.cbcb.ghost.align

import scopt._
import scala.io._
import scala.math._
import scala.collection.mutable._
import scala.util.Random
import scala.collection.JavaConversions._
import scala.tools.nsc.io._
import scala.concurrent.duration._
import scala.concurrent.Future

import akka.actor.{Actor, ActorRef, ActorSystem, Props}
import akka.pattern.ask
import akka.util.Timeout
// Akka 1.x
//import akka.actor.{ Actor, ActorRef, ActorRegistry }
//import akka.dispatch.{ Future, Dispatchers, MessageDispatcher }

import java.{ io => jio }
import java.io.{ FileWriter => JFileWriter, PrintStream => JPrintStream }

import org.jgrapht.graph._
import org.jgrapht.alg.{ ConnectivityInspector }
import org.ejml.simple._
import org.apache.mahout.math.{ DenseMatrix }
import org.apache.mahout.math.decomposer.lanczos.{ LanczosSolver }

import net.robpatro.utils.time.Timer
import net.robpatro.utils.console.ProgressBar
import edu.umd.cbcb.ghost.io.{ GEXFReader, EdgeListReader, SimpleVertex }
import edu.umd.cbcb.ghost.graphtools.actors._
import edu.umd.cbcb.ghost.graphtools.Utils
import edu.umd.cbcb.ghost.graphtools.{ LaplacianCalculator, SequenceWeightedEdge }

object ComputeSubgraphSignatures {

  type GraphT = WeightedPseudograph[SimpleVertex, SequenceWeightedEdge]
  type SubgraphT = UndirectedWeightedSubgraph[SimpleVertex, SequenceWeightedEdge]
  /**
   * Configuration Options
   */
  class Config {
    var numActors = 0
    var infile = ""
    var ghtree = ""
    var outdir = ""
    var numHops = 4
    var simpleSig = false
    var useWeight = false
  }

  def main(args: Array[String]) {
    var config = new Config
    val parser = new OptionParser("ComputeSubgraphSignatures") {
      intOpt("p", "numProcessors", "number of actors to use in parallel",
	   { v: Int => config.numActors = v })
      opt("i", "input", "input network file",
	{ v: String => config.infile = v})
      opt("g", "ghtree", "Gomory-Hu cut tree",
	{ v: String => config.ghtree= v})
      opt("o", "output", "output signature file",
	{ v: String => config.outdir = v})
      intOpt("k", "numHops", "maximum hop neighborhood to consider",
	     { v: Int => config.numHops = v})
      booleanOpt("s", "simpleSig", "output (simple) signature format",
		 {  v: Boolean => config.simpleSig = v} )
      booleanOpt("w", "useWeight", "use edge confidence",
	       { v: Boolean => config.useWeight = v})
    }

    if (parser.parse(args)) {

      val wsprefix = if (config.useWeight) { "" } else { "Not " }
      println(wsprefix+"using edge confidence")
      println(config.useWeight)

      // Resolve issue with Java7 and Scala 2.9.1 (SI-4960)
      val version = System.getProperty("java.runtime.version")
      if ( version.slice(0,3) == "1.7" ) { 
        System.setProperty("java.vm.vendor", "Sun Microsystems Inc.")
      }

      // Create a new EdgeListReader and read in the graph
      val elr = new GEXFReader(config.infile, config.useWeight)
      val H = elr.parse

      /* Build the Gomory-Hu Tree from file */
      val nameToVertexMap = H.vertexSet.map{ x => (x.name, x) }.toMap
      val gomoryHuTree = new SimpleWeightedGraph[SimpleVertex, SequenceWeightedEdge]( classOf[SequenceWeightedEdge] )
      if ( config.ghtree != "" ) {
	H.vertexSet.foreach{ v => gomoryHuTree.addVertex(v) }
	Source.fromFile( config.ghtree ).getLines.foreach{
	  l =>
	    val Array(s, t, c) = l.split(" ")
	    val e = new SequenceWeightedEdge
	    gomoryHuTree.addEdge( nameToVertexMap(s), nameToVertexMap(t), e )
	    gomoryHuTree.setEdgeWeight(e, c.toDouble )
	}
      }

      val maxFlowMat = if (config.ghtree != "") {  Some( Utils.gomoryHuTreeToMatrix( gomoryHuTree ) ) } else { None }

      // val ci = new ConnectivityInspector(G)Source
      // Sort the connected components in *descending* order
      // val vertexSets = ci.connectedSets.sortWith { (x, y) => x.size > y.size }

      // Grab the largest connected component of G
      // val H = new SubgraphT(G, vertexSets(0))

      import java.io.{ File => JFile, FileWriter => JFileWriter }
      import java.io.{ FileOutputStream => JFileOutputStream, FileInputStream => JFileInputStream }
      import java.io.{ BufferedOutputStream => JBufferedOutputStream }
      import java.io.{ DataOutputStream => JDataOutputStream }

      // Relabel the verticies
      H.vertexSet.zipWithIndex.foreach { vi => vi._1.id = vi._2 }
      println("\tLargest connected component, H : size = " + H.edgeSet.size + ", order = " + H.vertexSet.size)
      val reqNumSigs = H.vertexSet.size
      val k = config.numHops

      // Akka 1.x
      // Create a message dispatcher to load-balance our actors
      //val workStealingDispatcher = Dispatchers.newExecutorBasedEventDrivenWorkStealingDispatcher("subgraph-extraction-dispatcher")
      //.withNewThreadPoolWithLinkedBlockingQueueWithUnboundedCapacity
      //.setCorePoolSize(config.numActors)
      //.build

      val extractionActors = ArrayBuffer.empty[ActorRef]
      val system = ActorSystem("SubgraphExtractionActorSystem")
      // Create the actors that will be responsible for extracting the subgraphs
      def createSubgraphExtractors( rec : ActorRef, nex : Int ) = {
	(0 until nex).foreach{ i =>
          // Akka 1.x
	  //val actor = Actor.actorOf( new SubgraphExtractor(H, k, rec, maxFlowMat, workStealingDispatcher) )
          //actor.start
          // Akka 2.0
          val actor = system.actorOf( Props(new SubgraphExtractor(H, k, rec, maxFlowMat)), name="subgraph-extractor-%d".format(i) )
          extractionActors += actor
	}
      }

      // Needed to output the descriptors to file
      import java.io.{ File => JFile, FileWriter => JFileWriter }
      import java.io.{ FileOutputStream => JFileOutputStream, FileInputStream => JFileInputStream }
      import java.io.{ BufferedOutputStream => JBufferedOutputStream }
      import java.io.{ DataOutputStream => JDataOutputStream }
      import java.util.zip.{ GZIPOutputStream => JGZIPOutputStream, CheckedOutputStream => JCheckedOutputStream }

      var descFile : Either[ JDataOutputStream, JPrintStream ] = null.asInstanceOf[ Either[JDataOutputStream, JPrintStream] ]
      var gzOutStream : JGZIPOutputStream = null.asInstanceOf[ JGZIPOutputStream ]
      if (config.simpleSig) {
	descFile = new Right( new JPrintStream( config.outdir ) )
      } else {
	val dos = new JDataOutputStream( new JGZIPOutputStream ( new JFileOutputStream( config.outdir ) ) )
	descFile = new Left( dos )
      }

      // The actor the receives the subgraph descriptors and writes them to file
      // Akka 1.x
      //val descriptorReceiver = Actor.actorOf( new SubgraphDescriptorReceiver( reqNumSigs, k, descFile ) )
      //descriptorReceiver.start
      // Akka 2.x
      val descriptorReceiver = system.actorOf( Props(new SubgraphDescriptorReceiver( reqNumSigs, k, descFile )) )
      // Create the actual subgraph extraction actors
      createSubgraphExtractors( descriptorReceiver, config.numActors )
      val numExtractors = config.numActors

      // Populate their work queues randomly (the acutal load balancing will be done using the dispatcher created above)
      val rgen = new Random()
      // Akka 1.x
      //val extractionActors = Actor.registry.actorsFor[ SubgraphExtractor[GraphT] ]( classOf[SubgraphExtractor[GraphT]] )
      //val extractionActors = context.actorsFor[ SubgraphExtractor[GraphT] ]( classOf[SubgraphExtractor[GraphT]] )
      H.vertexSet.foreach{ v => extractionActors( rgen.nextInt(numExtractors) ) ! ExtractSubgraphDescriptor( v ) }

      // Until we've computed the necessary number of descriptors . . . wait
      while (!descriptorReceiver.isTerminated) {
        Thread.sleep(1000)
      }

      /*var dfuture: Future[Boolean]  = ask(descriptorReceiver, CompleteQuery())(2 seconds).mapTo[Boolean]
      while ( !( dfuture.isCompleted && dfuture.value.get.right.get ) ) {
        dfuture = ask(descriptorReceiver, CompleteQuery())(2 seconds).mapTo[Boolean]
        Thread.sleep(1500)
      }
      */

      // Akka 1.x
      //var isDone  = descriptorReceiver.ask(CompleteQuery())(2 seconds)
      //while ( !( isDone.isDefined && isDone.get) ) {
      //isDone = descriptorReceiver.ask(CompleteQuery())(2 seconds)
      //Thread.sleep(1500)
      //}

      // Once we're done, stop the descriptor receiver and all of the subgraph extractors
      descriptorReceiver ! Stop()
      system.stop(descriptorReceiver)
      system.shutdown()
      // Akka 1.x
      //Actor.registry.shutdownAll

      println(config.simpleSig)
      // Close the output file
      if (config.simpleSig) {
	descFile.right.get.close
      } else {
	descFile.left.get.close
      }

      System.exit(0)
    }

  }
}
