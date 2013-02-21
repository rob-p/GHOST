package edu.umd.cbcb.ghost.align

import scopt._

import scala.io._
import scala.math._
import scala.collection.mutable._
import scala.util.Random
import scala.collection.JavaConversions._
import scala.tools.nsc.io._

import java.{ io => jio }
import java.io.{ FileWriter => JFileWriter }

import org.jgrapht.graph._
import org.jgrapht.alg.{ ConnectivityInspector }
import org.apache.mahout.math.{ DenseMatrix }
import org.apache.mahout.math.decomposer.lanczos.{ LanczosSolver }

import net.robpatro.utils.time.Timer
import net.robpatro.utils.console.ProgressBar
import edu.umd.cbcb.ghost.io.{ GEXFReader, EdgeListReader, SimpleVertex }
import edu.umd.cbcb.ghost.graphtools.{ LaplacianCalculator, SequenceWeightedEdge }

object SpectralDecomposition {

  /**
   * Configuration Options
   */
  class Config {
    var numActors = 0
    var infile = ""
    var outdir = ""
    var seqfile = ""
    var numHops = 4
    var useWeight = false
  }


  def main(args: Array[String]) {
    var config = new Config
    val parser = new OptionParser("ComputeSubgraphSignatures") {
      intOpt("p", "numProcessors", "number of actors to use in parallel",
	   { v: Int => config.numActors = v })
      opt("i", "input", "input network file",
	{ v: String => config.infile = v})
      opt("s", "seqfile", "sequence input file",
	{ v: String => config.seqfile = v})
      opt("o", "output", "output signature file",
	{ v: String => config.outdir = v})
      intOpt("k", "numHops", "maximum hop neighborhood to consider",
	     { v: Int => config.numHops = v})
      booleanOpt("w", "useWeight", "use edge confidence",
	       { v: Boolean => config.useWeight = v})
    }

    if (parser.parse(args)) {
     val wsprefix = if (config.useWeight) { "" } else { "Not " }
      println(wsprefix+"using edge confidence")
      println(config.useWeight)

      // Create a new EdgeListReader and read in the graph
      val elr = new GEXFReader(config.infile, config.useWeight)
      val G = elr.parse

      val ci = new ConnectivityInspector(G)
      // Sort the connected components in *descending* order
      val vertexSets = ci.connectedSets.sortWith { (x, y) => x.size > y.size }

      // Grab the largest connected component of G
      val H = new Subgraph[SimpleVertex, SequenceWeightedEdge, WeightedPseudograph[SimpleVertex, SequenceWeightedEdge]](G, vertexSets(0))

      import java.io.{ File => JFile, FileWriter => JFileWriter }
      import java.io.{ FileOutputStream => JFileOutputStream, FileInputStream => JFileInputStream }
      import java.io.{ BufferedOutputStream => JBufferedOutputStream }
      import java.io.{ DataOutputStream => JDataOutputStream }

      println("SEQFILE: "+config.seqfile)
      if ( !config.seqfile.isEmpty ) {
	import org.biojava3.alignment.SmithWaterman._;
	import org.biojava3.alignment._;
	import org.biojava3.alignment.template.SequencePair;
	import org.biojava3.alignment.template.SubstitutionMatrix;
	import org.biojava3.core.sequence.ProteinSequence;
	import org.biojava3.core.sequence.compound.AminoAcidCompound;
	import org.biojava3.core.sequence.io.FastaReaderHelper;

	val pb = new ProgressBar( H.edgeSet.size, "=" )

	val seqMap = FastaReaderHelper.readFastaProteinSequence( new JFileInputStream(config.seqfile) )

	val subMat = new SimpleSubstitutionMatrix[AminoAcidCompound]
	val aligner = new SmithWaterman[ProteinSequence, AminoAcidCompound]
	aligner.setGapPenalty( new SimpleGapPenalty )
	aligner.setSubstitutionMatrix( subMat )
	var ctr = 0
	H.edgeSet.foreach{
	  eij =>
	    val (srcKey, tgtKey) = ( G.getEdgeSource(eij).name, G.getEdgeTarget(eij).name )
            val (s0, s1) = ( seqMap(srcKey), seqMap(tgtKey) )
	    aligner.setQuery(s0); aligner.setTarget(s1)
	    val sim = aligner.getSimilarity
	    eij.seqWeight = sim
	    pb.update(ctr); ctr += 1
	}

	pb.done

      }

      // Relabel the verticies
      H.vertexSet.zipWithIndex.foreach { vi => vi._1.id = vi._2 }
      println("\tLargest connected component, H : size = " + H.edgeSet.size + ", order = " + H.vertexSet.size)

      val timer = new Timer; timer.start

      val lc = new LaplacianCalculator(H, (x, y) => y )
      val (lap, vols) = lc.compute

      timer.stop
      println("Building Laplacian took "+timer.elapsed / 1000.0+" seconds")
      timer.reset

      val N = lap.numRows
      val M = N//min(300, N-1)

      /** Use the sparse eigensolver **/
      //val solver = new LanczosSolver
      //var evals = new java.util.ArrayList[java.lang.Double]
      //val eigens = new DenseMatrix(M, N)
      //solver.solve(lap, M, eigens, evals, true)

      /** Use the dense eigensolver **/
      import no.uib.cipr.matrix._
      val denseL = new UpperSPDDenseMatrix(N)
      for (i <- 0 until N;
	   j <- 0 until N; if j >= i) {
	denseL.set(i, j, lap.get(i,j))
	denseL.set(j, i, lap.get(j,i))
      }


      val matFileName = "mat.bin"
      var tstream = new JFileOutputStream(matFileName)
      var tbuffer = new JDataOutputStream(new JBufferedOutputStream(tstream))
      // Write the matrix size
      tbuffer.writeInt(N); tbuffer.writeInt(N)
      // Write the entries
      denseL.iterator.foreach {  x => tbuffer.writeDouble(x.get()) }
      tbuffer.close()


      val sdEVD = new SymmDenseEVD(N, true)
      println("Solving Eigenproblem")
      timer.start
      sdEVD.factor(denseL)
      val devals = sdEVD.getEigenvalues()
      val devecs = sdEVD.getEigenvectors()
      timer.stop
      println("Solving the Eigenproblem for "+M+" pairs took "+timer.elapsed/1000.0+" seconds ")
      devals.zipWithIndex.foreach{ xi => print("evals["+xi._2+"] = "+xi._1+"\n")}

      val outDir = new Directory(new java.io.File(config.outdir))
      outDir.createDirectory(false, false)


      val specFileName = List(outDir, "spectrum.txt").mkString(File.separator)
      //val specFile = new JFile(specFileName)
      var stream = new JFileOutputStream(specFileName)
      var buffer = new JDataOutputStream(new JBufferedOutputStream(stream))
      // This conversion is *UGLY*, is there another way to do it?
      devals.iterator.foreach { x => buffer.writeDouble(x.asInstanceOf[Double]) }
      buffer.close()


      val evecFileName = List(outDir, "evecs.txt").mkString(File.separator)
      stream = new JFileOutputStream(evecFileName)
      buffer = new JDataOutputStream(new JBufferedOutputStream(stream))
      // Write the matrix size
      buffer.writeInt(N); buffer.writeInt(M)
      // Write the entries
      var d = 0.0
      devecs.iterator.foreach {  x => buffer.writeDouble(x.get()) }
      buffer.close()

      val volFileName = List(outDir, "vols.txt").mkString(File.separator)
      stream = new JFileOutputStream( volFileName )
      buffer = new JDataOutputStream( new JBufferedOutputStream(stream) )
      // Write the entries
      vols.iterator.foreach {  x => buffer.writeDouble(x.get()) }
      buffer.close()

      val idxMapFileName = List(outDir, "idxMap.txt").mkString(File.separator)
      val idxMapFile = new JFileWriter( idxMapFileName )
      H.vertexSet.foreach{
	vi =>
	  idxMapFile.write(vi.id+"\t"+vi.name+"\n")
      }
      idxMapFile.close

      /** Print the Laplacian
      for( i <- 0 until lap.rowSize; j <- 0 until lap.columnSize ) {
       	val e = lap.get(i,j)
       	if (e != 0) {
       	  println("L["+i+","+j+"] = "+lap.get(i,j))
       	}
      }
      **/

      System.exit(1)
    }

  }
}
