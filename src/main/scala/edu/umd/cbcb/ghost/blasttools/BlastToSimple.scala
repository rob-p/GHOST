package edu.umd.cbcb.ghost.blasttools

import scala.collection.mutable.{ OpenHashMap => Map }
import scala.math._
import scalax.io._
import scalax.file._
import org.rogach.scallop._
import org.rogach.scallop.exceptions._

object BlastToSimple {

  def parseSelfBlast( fname: String, cidx: Int ) = {
    val m = Map.empty[String, Double]
    val ilines = Resource.fromFile( fname ).lines()
    ilines.foreach{ l =>
      val lt = l.trim
       if ( !lt.startsWith("#") ) {
         val toks = lt.split("""\s+""")
         if (toks(0) == toks(1)) { m(toks(0)) = toks(cidx).toDouble }
       }
    }
    m
  }

  def parsePairwiseBlast( fname: String, cidx: Int ) = {
    val m = Map.empty[(String, String), Double]
    val ilines = Resource.fromFile( fname ).lines()

    ilines.foreach{ l =>
      val lt = l.trim
      if (!lt.startsWith("#")) {
        val toks = lt.split("""\s+""")
        m( (toks(0), toks(1)) ) = toks(cidx).toDouble
      }
    }
    m
  }

  def main( argv: Array[String] ) {
    
    object Conf extends ScallopConf(argv) { //Seq("--help")) {
      version("BlastToSimple v1.0.0 (c) 2012 Rob Patro")
      banner("""Usage: BlastToSimple [OPTION]... [tree|palm] [OPTION]... [tree-name]
               |test is an awesome program, which does something funny      
               |Options:
               |""".stripMargin)
      footer("\nFor all other tricks, consult the documentation!")
      // ... options ...
      val lMap = opt[String]("lscores", descr="self-BLAST scores of left hand side", required = false)
      val rMap = opt[String]("rscores", descr="self-BLAST scores of right hand side", required = false)
      val psMap = opt[String]("psScores", descr="pairwise BLAST scores")
      val cidx = opt[Int]("cidx", descr="score index", required=true)
      val ofile = opt[String]("ofile", descr="output file", required = true )

      codependent( lMap, rMap )
      verify
    }

    var omap = Map.empty[ (String, String), Double ]
    var ofile = Option.empty[Path]

    val cidx = Conf.cidx()    
    assert( cidx > 1 && cidx < 12 )
    val psMap = parsePairwiseBlast(Conf.psMap(), cidx)    
    ofile = Some(Path.fromString(Conf.ofile()))
      
    if ( !Conf.lMap.get.isEmpty ) {
      val lMap = parseSelfBlast(Conf.lMap(), cidx)
      val rMap = parseSelfBlast(Conf.rMap(), cidx)
      psMap.foreach{
        case ((k, hk), d) => {
          try {
            omap((k,hk)) = d / ( sqrt(1.0 - min(lMap(k),0.999)) * sqrt(1.0 - min(0.999,rMap(hk))) )
          } catch {
            case e: Throwable => println(e)
          }
        }
      }
    } else {
      omap = psMap
    }
    
    // argv.size match { 
    //   case 3 => { 

    //   val cidx = argv(1).toInt
    //   assert( cidx > 1 && cidx < 12 )

    //   omap = parsePairwiseBlast( argv(0), cidx )
    //   ofile = Some(Path.fromString( argv(2) ))
    //   }
    //   case 5 => { 

    //   val cidx = argv(3).toInt
    //   assert( cidx > 1 && cidx < 12 )

    //   val lMap = parseSelfBlast( argv(0), cidx )
    //   val rMap = parseSelfBlast( argv(1), cidx )
    //   val psMap = parsePairwiseBlast( argv(2), cidx )
    //   ofile = Some(Path.fromString( argv(4) ))
    //   psMap.foreach{
    //     case ((k, hk), d) => {
    //       try {
    //         omap((k,hk)) = d / ( sqrt(1.0 - min(lMap(k),0.999)) * sqrt(1.0 - min(0.999,rMap(hk))) )
    //       } catch {
    //         case e => println(e)
    //       }
    //     }
    //   }
    //   }
    //   case _ => { 
    //     println("usage: BlastToSimple ifile scoreidx ofile\nBlastToSimple lscores rscores ifile scoreidx ofile")
    //   }
    // }

    ofile.get.createFile( failIfExists=false )
    for{
       // create a processor (signalling the start of a batch process)
       processor <- ofile.get.outputProcessor
       // create an output object from it
       out = processor.asOutput
    }{
        // all writes to out will be on the same open output stream/channel
        omap.foreach{
          case ((l,r), d) => {
            out.write(l+"\t"+r+"\t"+d+"\n")
        }
      }
    }

  }

}
