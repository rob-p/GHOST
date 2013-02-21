package edu.umd.cbcb.ghost.align

import akka.actor.ActorSystem
import com.typesafe.config.ConfigFactory

import scopt._
import scala.io._
import scala.math._
import scala.collection.mutable.{ HashSet => MHashSet, HashMap => MHashMap, Set => MSet, OpenHashMap => MMap, ListBuffer }
import scala.collection.parallel.{ ParIterable, ParMap => PMap }
import scala.collection.parallel.mutable.{ ParHashSet => PMSet, ParHashMap => PMHashMap, ParArray }
import scala.util.Random
import scala.collection.JavaConversions._
import scala.tools.nsc.io._
import scala.io.Source
import scala.collection.mutable.PriorityQueue

import java.{ io => jio }
import java.io.{ FileWriter => JFileWriter, BufferedWriter => JBufferedWriter, PrintWriter => JPrintWriter }
import java.util.{ HashMap => JHashMap }

import org.jgrapht.graph._
import org.jgrapht.alg.{ ConnectivityInspector }
import org.ejml.simple.SimpleMatrix
import org.ejml.ops.NormOps

import net.robpatro.utils.time.Timer
import net.robpatro.utils.console.ProgressBar
import edu.umd.cbcb.ghost.io.{ GEXFReader, EdgeListReader, SimpleVertex }
import edu.umd.cbcb.ghost.graphtools.{ Utils, SequenceWeightedEdge }
import edu.umd.cbcb.ghost.graphtools.actors._
import edu.umd.cbcb.ghost.align.alignmentextenders._
import edu.umd.cbcb.ghost.mathutils.TFIDFWeighter
import edu.umd.cbcb.ghost.align.matchers._
import edu.umd.cbcb.ghost.align.localsearch._
import edu.umd.cbcb.ghost.mathutils.UtilFunctions

import grizzled.config.Configuration

object LocalAligner{

  val logoString = """
     ______ __  __ ____  _____ ______
    / ____// / / // __ \/ ___//_  __/
   / / __ / /_/ // / / /\__ \  / /
  / /_/ // __  // /_/ /___/ / / /
  \____//_/ /_/ \____//____/ /_/
  """
  type VertexT = SimpleVertex
  type EdgeT = SequenceWeightedEdge
  type GraphT = WeightedPseudograph[SimpleVertex, SequenceWeightedEdge]
  type SubgraphT = Subgraph[SimpleVertex, SequenceWeightedEdge, GraphT]

  /**
   * Configuration Options
   */
  class Config {
    var numActors = 2
    var numBins = 20
    var localSearchIter = Option.empty[Int]//5
    var constrainedRatio = Option.empty[Double]//1.0
    var alpha = Option.empty[Double]
    var beta = Option.empty[Double]
    var g1FileName = ""
    var g2FileName = ""
    var g1SigFile = ""
    var g2SigFile = ""
    var blastFile = ""
    var outputDir = ""
    var cfgFile = ""
  }

  def CreateGraph(sg: SubgraphT) = {
    var vnmap = sg.vertexSet.map{ v => (v, (v.name,v.id)) }.toMap
    var nvmap = new MHashMap[String, VertexT]
    val g = new GraphT(classOf[ EdgeT ] );
    sg.vertexSet.foreach{
      v =>
	val nid = vnmap(v);
      val nvert = new VertexT(nid._1, nid._2);
      nvmap(nid._1) = nvert
      g.addVertex(nvert)
    };
    sg.edgeSet.foreach{
      e =>
	val src = nvmap( sg.getEdgeSource(e).name )
      val tgt = nvmap( sg.getEdgeTarget(e).name )
      val we = new SequenceWeightedEdge
      g.addEdge( src, tgt, we )
      g.setEdgeWeight( we, 1.0 )
    }
    g
  }

  def smallestKElems[CT](xs: TraversableOnce[CT], k: Int, sf: (CT) => Double ) = {
      var ss = scala.collection.SortedSet.empty[Double]
      var min = Double.NegativeInfinity
      var len = 0
      xs foreach { elem =>
        val e = sf(elem)
        if (ss.size < k || e < min) {
          ss = ss + e
          min = ss.last
         }
         if (ss.size > k) {
            ss = ss - ss.last
            min = ss.last
         }                    
      }
      ss
  } 

  def computeDistances(
    lSigMap: SignatureMap,
    rSigMap: SignatureMap,
    unmatchedLeft: ParIterable[Int],
    unmatchedRight: ParIterable[Int],
    blastScores: Option[ScoreMap],
    a: Double,
    k: Int = 5,
    inferAlpha: Boolean ) = {

    val haveBlastScores = blastScores.isDefined
    val alpha = if(haveBlastScores) { a } else { 1.0 }
    val seqDist = blastScores.getOrElse( ScoreMap(None, 0.0) )
    val khop = lSigMap.maxHop
    val (maxG, maxH) = (lSigMap.maxDeg, rSigMap.maxDeg)
    val pb = new ProgressBar( unmatchedLeft.size, "=" )
    var i = 0

    val matches = unmatchedLeft.flatMap{
      lind: Int =>
	val lname = lSigMap.idToName(lind)
        var lspec = lSigMap.spectra(lname)
        
      val dlist = unmatchedRight.map{
	rind: Int =>
	val rname = rSigMap.idToName(rind)
        var rdeg = rSigMap.G.degreeOf( rSigMap.nameVertexMap(rname) )
        var rspec = rSigMap.spectra(rname)

        val structDist = UtilFunctions.structDistance( lspec, rspec )
        var dist = (structDist, seqDist((lname,rname)))

	(lname, rname, dist)
      }

      pb.update(i); i += 1
      dlist
    }

    pb.done()

    /*
    import java.io._
    {
      val obuf = new PrintWriter(new File("distancesStruct.txt" ))
      matches.seq.foreach{ case (l,r,(std, seqd)) => obuf.write(std+" ") }
      obuf.close()
    }

    {
      val obuf = new PrintWriter(new File("distancesSeq.txt" ))
      matches.seq.foreach{ case (l,r,(std, seqd)) => obuf.write(seqd+" ") }
      obuf.close()
    }
    */
   
    type matchType = (String, String, (Double, Double))
    def sortByStruct( m: matchType ) = m match { case (_,_,(d,_)) => d }
    def structCmp( m1: matchType, m2: matchType ) = { m1._3._1 < m2._3._1 }
    def sortBySeq( m: matchType ) = m match { case (_,_,(_,d)) => d }
    def seqCmp( m1: matchType, m2: matchType ) = { m1._3._2 < m2._3._2 }
    
    println("Evaluated %d potential matches".format(matches.size))
    if (inferAlpha) {


      val idx = (matches.seq.toSeq.length * 0.0015).toInt
      println("index = %d".format(idx))

      val sortedStructs = matches.seq.toSeq.sortWith( structCmp ).filter( x => x._3._1 > 0.0 )
      val stIdx = (sortedStructs.size * 0.0005).toInt
      val structVal = sortedStructs(stIdx)._3._1

      val sortedSeqs = matches.seq.toSeq.sortWith( seqCmp ).filter( x => x._3._2 > 0.0 )
      val seIdx = (sortedSeqs.size * 0.0005).toInt
      val seqVal = sortedSeqs(seIdx)._3._2

      println(s"Structural Value = $structVal, Sequence Value = $seqVal")

      val cAlpha = seqVal / structVal
      println(s"topStruct = $structVal, topSeq = $seqVal, cAlpha = $cAlpha")
      
      //val cAlpha = 5e-12
      (matches, cAlpha)
    } else {
      (matches, alpha)
    }
    

  }

  def getNearestNeighbors(
    matches: ParIterable[(String, String, (Double, Double))],
    alpha: Double ) = {

    /* Create the priority queue; ensuring the implicit ordering above
     * is used.  Insert the scaled matches into the queue
     */
    val pqord = Ordering[Double].on[(String, String, Double)]{ m => 1.0 - m._3 }

    // Create a map (PID1, PID2) => Score
    def combineWithAlpha( dists: (Double, Double) ) = {
      alpha * dists._1 + (1.0-alpha) * dists._2
    }

    val distMatches = matches.map{ m => (m._1, m._2, combineWithAlpha(m._3) ) }.toList
    Random.setSeed(System.currentTimeMillis)
    Random.shuffle(distMatches)
    val matchMap = MHashMap.empty[(String,String), Double]
    distMatches.foreach{ case (l,r,d) => matchMap((l,r)) = d }

    // As well a priority Queue
    val mset = PriorityQueue()(pqord)
    mset ++= distMatches
    println("priority of of pq.head should be close to 0.0; it is "+mset.head._3)
    println("There are %d potential matches in the queue".format(mset.size))
    (matchMap, mset)
  }
    
  
  def computeNearestNeighbors(
    lSigMap: SignatureMap,
    rSigMap: SignatureMap,
    unmatchedLeft: ParIterable[Int],
    unmatchedRight: ParIterable[Int],
    blastScores: Option[ScoreMap],
    a: Double,
    k: Int = 5 ) = {

    val haveBlastScores = blastScores.isDefined
    val alpha = if(haveBlastScores) { a } else { 1.0 }
    val seqDist = blastScores.getOrElse( ScoreMap(None, 0.0) )
    val khop = lSigMap.maxHop
    val (maxG, maxH) = (lSigMap.maxDeg, rSigMap.maxDeg)

    val pb = new ProgressBar( unmatchedLeft.size, "=" )
    var i = 0

    val matches = unmatchedLeft.flatMap{
      lind: Int =>
	val lname = lSigMap.idToName(lind)
        var lspec = lSigMap.spectra(lname)
        
      val dlist = unmatchedRight.map{
	rind: Int =>
	val rname = rSigMap.idToName(rind)
        var rdeg = rSigMap.G.degreeOf( rSigMap.nameVertexMap(rname) )
        var rspec = rSigMap.spectra(rname)

        val structDist = UtilFunctions.structDistance( lspec, rspec )

        var dist = alpha * structDist

        // Feb 20 -- histogram dist
	//var dist = alpha * (invFeatWeight * UtilFunctions.l1Distance( wvec, lvec, rvec ) )

        // !!!! Sequence Distance
	dist += (1.0 - alpha) * seqDist( (lname, rname) )

        //} else {
        //  dist += (1.0 - alpha)
        //}
	// val density = 0.5 * (ldeg + rdeg) / (lSigMap.maxDeg + rSigMap.maxDeg).toFloat
        // dist *= (1.0 - density)
	(rname, dist)
      }.toList
      /*
      val numMatches = dlist.size
      val dlist2 = dlist.sortBy{ x => x._2 }.view.zipWithIndex.map{ case ((rname, dist), ind) => 
        val relTopRank = ind.toFloat / numMatches
        val relSeqRank = seqDist( (lname, rname) )
        (rname, alpha*relTopRank + (1.0-alpha)*relSeqRank)
      }.toList
      */
      pb.update(i); i += 1

      val dists = Random.shuffle(dlist).toArray

      val ndists = dists.sortBy{ x => x._2 }.slice(0,k).map{ m => (lname, m._1, m._2) }
      ndists
    }.seq

    pb.done()

    println("There are %d potential matches here".format(matches.size))

    /* Create the priority queue; ensuring the implicit ordering above
     * is used.  Insert the scaled matches into the queue
     */
    val pqord = Ordering[Double].on[(String, String, Double)]{ m => 1.0 - m._3 }

    Random.setSeed(System.currentTimeMillis)
    Random.shuffle(matches)

    //System.err.println( matches.map{ x => x._3 }.mkString(" ") )
    // Create a map (PID1, PID2) => Score
    val matchMap = matches.map{ m => ( (m._1, m._2), m._3 ) }.toMap
    // As well a priority Queue
    val mset = PriorityQueue()(pqord)
    mset ++= matches
    println("priority of of pq.head should be close to 0.0; it is "+mset.head._3)
    println("There are %d potential matches in the queue".format(mset.size))
    (matchMap, mset)
  }

  def main(args: Array[String]) {
    var config = new Config
    val parser = new OptionParser("ComputeSubgraphSignatures") {
      intOpt("p", "numProcessors", "number of actors to use in parallel",
	     { v: Int => config.numActors = v})
      intOpt("n", "numBins", "number of histogram bins to use",
             { v: Int => config.numBins = v })
      intOpt("l", "localSearch", "number of local search iterations to perform",
             { v: Int => config.localSearchIter = Some(v) })
      doubleOpt("a", "alpha", "alpha parameter (trade-off between sequence & topology)",
             { v: Double => config.alpha = Some(v) })
      doubleOpt("b", "beta", "beta parameter (Blast e-value cutoff)",
             { v: Double => config.beta = Some(v) })
      doubleOpt("r", "unconstrained", "percentage of local swaps which should be allowed to be unconstrained",
             { v: Double => config.constrainedRatio = Some(v) } )
      opt("u", "g1", "first input graph",
	{ v: String => config.g1FileName = v.trim})
      opt("v", "g2", "second input graph",
	{ v: String => config.g2FileName = v.trim})
      opt("s", "s1", "first input signature file",
	{ v: String => config.g1SigFile = v.trim})
      opt("t", "s2", "second input signature file",
	{ v: String => config.g2SigFile = v.trim})
      opt("b", "blast", "Pairwise BLAST Scores",
	{ v: String => config.blastFile = v.trim})
      opt("o", "output", "output signature file",
	{ v: String => config.outputDir = v.trim})
      opt("c", "config", "configuration file",
	{ v: String => config.cfgFile = v.trim})
      help("h", "help", "produce this help message")
    }

    if (parser.parse(args)) {

      // Resolve issue with Java7 and Scala 2.9.1 (SI-4960)
      val version = System.getProperty("java.runtime.version")
      if ( version.slice(0,3) == "1.7" ) { 
        System.setProperty("java.vm.vendor", "Sun Microsystems Inc.")
      }

      
      //scala.concurrent.context.numThreads = config.numActors
      
      val cfg = Configuration( Source.fromFile(config.cfgFile) )
      val mainCfgOpt = cfg.getSection("main")
      assert( !mainCfgOpt.isEmpty, "Configuration file must contain [main] section" )
      val mainCfg = mainCfgOpt.get

      println(logoString)

      //println("Options :")
      //mainCfg.options.foreach{ kv => println("\t %s : %s".format(kv._1, kv._2)) }

      def getOrDie(k: String, msg: String) = {
	val v = mainCfg.options.get(k)
	assert(v.isDefined, msg)
	v.get
      }
      /* Required configuration variables */
      val g1File = getOrDie("network1", "configuration must have key \"network1\"")
      val g2File = getOrDie("network2", "configuration must have key \"network2\"")

      val g1SigFile = getOrDie("sigs1", "configuration must have key \"sigs1\"")
      val g2SigFile = getOrDie("sigs2", "configuration must have key \"sigs2\"")

      /* Optional configuration variables */
      val numBinsStr = mainCfg.options.get("numbins").getOrElse("15"); val numBins = numBinsStr.toInt
      val blastFileOpt = mainCfg.options.get("sequencescores")
      val matchType = mainCfg.options.get("matcher").getOrElse("linear").trim()
      val beta = config.beta.getOrElse(mainCfg.options.get("beta").getOrElse("Infinity").toDouble)
      val constrainedRatio = config.constrainedRatio.getOrElse(mainCfg.options.get("ratio").getOrElse("1.0").toDouble)
      val localSearchIter = config.localSearchIter.getOrElse(mainCfg.options.get("searchiter").getOrElse("10").toInt)
      //val constrainedRatio = config.constrainedRatio
      val validMatchTypes = Set("linear", "quadratic")
      assert( validMatchTypes contains matchType, "matcher must be one of {linear | quadratic}")

      val confFile = new java.io.File("conf/ghost.conf")
      val customConf = ConfigFactory.parseFile( confFile )
      val system = ActorSystem("QAPExtractionActors", customConf)

      val matcherCreator = matchType match {
	case "linear" => Some(LAPMatcherCreator)
	case "quadratic" => Some(SpectralMatcherCreator)
	case _ => None
      }

      val nn = mainCfg.options.get("nneighbors")
      // Create a new GEXFReader and read in the graphs
      val g1f = new GEXFReader(g1File, false).parse
      val g2f = new GEXFReader(g2File, false).parse

      var sg1 = new SubgraphT(g1f, g1f.vertexSet) 
      var sg2 = new SubgraphT(g2f, g2f.vertexSet) 

      val g1 = CreateGraph(sg1)
      val g2 = CreateGraph(sg2)
      val (n1, n2) = (g1.vertexSet.size, g2.vertexSet.size)

      val avgDeg1 = g1.vertexSet.toList.map{ x => g1.degreeOf(x) }.sorted.apply(g1.vertexSet.size / 2 )
      val avgDeg2 = g2.vertexSet.toList.map{ x => g2.degreeOf(x) }.sorted.apply(g2.vertexSet.size / 2 )
      val avgDeg = 0.5 * ( avgDeg1 + avgDeg2 )
      

      val (sigMap1, sigMap2) = (new SignatureMap( g1SigFile, g1, numBins ), new SignatureMap( g2SigFile, g2, numBins ));

      val histBase = avgDeg
      print("Loading "+g1SigFile+" and computing signatures . . . ")
      sigMap1.readSignatures( g1SigFile, histBase )
      println("Done")
      print("Loading "+g2SigFile+" and computing signatures . . . ")
      sigMap2.readSignatures( g2SigFile, histBase )
      println("Done")

      var blastScores: Option[ ScoreMap ] = None
      var (maxScore, minScore) = (-Double.MaxValue, Double.MaxValue)
      if ( blastFileOpt.isDefined ) {
	println( "Using sequence information\nReading BLAST Scores from %s".format(blastFileOpt.get) )
        val maxSeqScore = beta
        val lit = Source.fromFile(blastFileOpt.get).getLines()
        val smap = collection.mutable.OpenHashMap.empty[(String, String), Double]
        val seqs = collection.mutable.ArrayBuffer.empty[(String, String, Double)]
	lit.foreach{
	  l =>
	    // get the input
	    val Array(s,t,v) = l.split('\t')
	    // keep track of the highest score
	    val d = min(maxSeqScore, v.toDouble)
	    minScore = min(minScore, v.toDouble)
            maxScore = max(maxScore, d)
	    // yield the map values
            if ( d < maxSeqScore ) {
              smap( (s,t) ) = d
              //seqs += ((s,t,d))
            }
        }
        println(s"MaxScore $maxScore, MinScore $minScore")
        /*
        // Insert sequences into the map by relative rank
        val seqGroups = seqs.groupBy{ x => x._1 }
        seqGroups.foreach{ case (prot1, matches) =>
          val numMatches = matches.size
          matches.sortBy{ x => x._3 }.view.zipWithIndex.foreach{ 
            case (m, i) => 
              smap( (prot1, m._2) ) =  i.toFloat/numMatches
          }
        }
        */
        //smap.foreach{ case (k,v) => smap(k) = v / maxScore } 
        smap.foreach{ case (k,v) => smap(k) = v / maxSeqScore }
	blastScores = Some( ScoreMap(Some(smap.par), if(beta == Double.PositiveInfinity){ 1.0 } else { 1.1 } ) )
      } else {
	println( "Not using sequence information" )
      }

      config.alpha = config.alpha match { 
        // If there was no alpha given on the cmdline, check for
        // one in the config file
        case None => mainCfg.options.get("alpha") match { 
          case Some(a) => Some(a.toDouble)
          case _ => Option.empty[Double]
        }
        case x => x
      }

      var (alpha, inferAlpha) = config.alpha match { 
        case None => { 
          blastScores match { 
            // If there are no BLAST scores
            case None => { (1.0,false) }
            // If there are BLAST scores
            case Some(scoreMap) => { (0.5,true) }
          }
        }
        // If the user provided an alpha explicitlly
        case Some(a) => (a,false)
      }
      println("inferAlpha = "+inferAlpha)
      //println("alpha = "+alpha)


      // Tried using mutable parallel hash set, but that collection
      // is currently broken (see SI-4678)
      var unmatchedLeft = Range(0, sigMap1.size).toSet.par
      var unmatchedRight = Range(0, sigMap2.size).toSet.par

      var matchedLeft = MHashSet.empty[String]
      var matchedRight = MHashSet.empty[String]

      var k = nn match {
	case Some("all") => sigMap1.size // unmatchedLeft.size
	case Some(n) => n.toInt
        case _ => 5
      }

      println( "Considering %d nearest neighbors".format(k) )
      
      val alignment = Alignment.empty[SimpleVertex]
      val currentAlignment = MMap.empty[SimpleVertex, (SimpleVertex,Double)]
      val fmap = MMap.empty[SimpleVertex, SimpleVertex]
      val rmap = MMap.empty[SimpleVertex, SimpleVertex]

      val pbar = new ProgressBar( sigMap1.size, "#" )
      val matches = MMap.empty[String, String]//MHashSet.empty[(String, String, Double)]
      var round = 0
      val (dmatches, cAlpha) = computeDistances(sigMap1, sigMap2, unmatchedLeft, unmatchedRight, blastScores, alpha, k, inferAlpha)
      var (mm, pq) = getNearestNeighbors(dmatches, alpha)//computeNearestNeighbors(sigMap1, sigMap2, unmatchedLeft, unmatchedRight, blastScores, alpha, k)

      while ( !(pq.isEmpty) && alignment.size < n1 ) { 
        val tmatch = pq.dequeue
	val (lm, rm, d) = tmatch
        val (lmv, rmv) = (sigMap1.nameVertexMap(lm), sigMap2.nameVertexMap(rm))

	if ( !alignment.isEitherAligned(lmv,rmv) ) { 
	  val numUnmatched = n1 - alignment.size
	  k = min(k, numUnmatched)
          
	  val extender = new AlignmentExtender(tmatch, sigMap1, sigMap2, blastScores,
                                               alpha, matchType, k, mm, alignment, pbar, matcherCreator, system)

          val lmatches = extender.extendAlignment
	}
      }

      // If we still have not aligned everything, do it radomly
      if ( alignment.size < n1 ) {
        alpha = 1.0
        inferAlpha = false
        val remLeft = sigMap1.idToName.collect{ case (id:Int,name:String) if(!alignment.isAlignedLeft(sigMap1.nameVertexMap(name))) => id }.toSet
        var remRight = sigMap2.idToName.collect{ case (id:Int,name:String) if(!alignment.isAlignedRight(sigMap2.nameVertexMap(name))) => id }.toSet

        remLeft.foreach{ x =>
          val y = remRight.head 
          alignment.insert( sigMap1.nameVertexMap( sigMap1.getName(x)) , sigMap2.nameVertexMap(sigMap2.getName(y)),1.0)
          remRight -= y
        }
        /*
        val (dmatches, cAlpha) = computeDistances(sigMap1, sigMap2, unmatchedLeft, unmatchedRight, None, alpha, k, inferAlpha)
        var (mm, pq) = getNearestNeighbors(dmatches, alpha)//computeNearestNeighbors(sigMap1, sigMap2, unmatchedLeft, unmatchedRight, blastScores, alpha, k)
              while ( !(pq.isEmpty) && alignment.size < n1 ) { 
                val tmatch = pq.dequeue
	        val (lm, rm, d) = tmatch
                val (lmv, rmv) = (sigMap1.nameVertexMap(lm), sigMap2.nameVertexMap(rm))

	        if ( !alignment.isEitherAligned(lmv,rmv) ) { 
	          val numUnmatched = n1 - alignment.size
	          k = min(k, numUnmatched)
            
	          val extender = new AlignmentExtender(tmatch, sigMap1, sigMap2, blastScores,
                                     alpha, matchType, k, mm, alignment, pbar, matcherCreator, system)

                  val lmatches = extender.extendAlignment
	       }
             }  
        */
      }
      
      alignment.foreach{ 
        case (lk, rk, d) =>
	fmap(lk) = rk
	rmap(rk) = lk
      }

      val finalMatches = alignment.map{ case (lk,rk,d) => (lk,rk) }
      val numMatches = finalMatches.filter{ m => m._1.name == m._2.name }.size
      val numPotentialMatches = finalMatches.size.toFloat // g1.vertexSet.size.toFloat
      println("Node Accuracy is "+numMatches+"/"+numPotentialMatches+" = "+100.0*(numMatches/numPotentialMatches)+"%")

      val gmap = new DefaultGraphMapping(fmap, rmap, g1, g2)
      val mappedVerts = finalMatches.map{ m => m._1 }.toSet
      var mappedSG = new SubgraphT(g1, mappedVerts)
      val e1 = mappedSG.edgeSet
      val nedges = e1.size
      val oedges = g2.edgeSet.size
      val nmapped = e1.count{ e1 => gmap.getEdgeCorrespondence(e1,true) != null }

      /*val inducedSubgraph = new GraphT( classOf[EdgeT] )
      val nameVertMap = MMap.empty[String, SimpleVertex]
      sg2.vertexSet.foreach{ v => inducedSubgraph.addVertex(v); nameVertMap(v.name) = v }
      inducedEdges.foreach{ e => val ne = new EdgeT(); inducedSubgraph.addEdge( nameVertMap(e.u), nameVertMap(e.v), ne ) } 
      */

      //println("ICS = %f %%".format(Utils.inducedConservedStructure(alignment, g1, g2)))

      val neb = nmapped
      //println("Edge Accuracy is %d / %d = %f %%".format(nmapped, nedges, 100.0*(nmapped.toDouble/nedges)) )
      //println("ICS = %f %%".format(Utils.inducedConservedStructure(alignment->-!, g1, g2)))
      //println("Phylo Dist is %d / %d = %f %%".format(nmapped, nedges, 1.0 - (nmapped.toDouble/min(nedges,oedges)) ))
      val ls = new NaiveSearch(alignment, g1, g2, blastScores)
      val nonConstrainedRatios = (0 until localSearchIter).map{ x => scala.math.exp(-x) }
      val totalSum = nonConstrainedRatios.sum 
      val ncInvSum = if (totalSum > 0.0) { constrainedRatio / totalSum } else { Double.PositiveInfinity }

      (0 until localSearchIter).foreach{ i => 
        val ncRatio = nonConstrainedRatios(i) * ncInvSum 
        println("%f %% of the search will be unconstrained".format(ncRatio))
        ls.searchIter( ncRatio )
      }

      //println("Edge Accuracy was %d / %d = %f %%".format(neb, nedges, 100.0*(neb.toDouble/nedges)) )
      val nea = Utils.mappedEdges( alignment.->-!, g1, g2 ).size
      println(s"Edge Accuracy $nea / $nedges = ${100.0*(nea.toDouble/nedges)}%")
      println(s"ICS = ${Utils.inducedConservedStructure(alignment->-!, g1, g2) * 100.0}%")
      // Write the alignment to file
      if (! config.outputDir.isEmpty ) {
       val ofile = new JPrintWriter( new JBufferedWriter( new JFileWriter(config.outputDir)) )
       alignment.foreach{ case (l,r,s) => ofile.println("%s\t%s".format(l.name, r.name)) }
       ofile.close
      }
      // Shutdown the actor system cleanly
      system.shutdown
    }
    System.exit(0)
  }

}

