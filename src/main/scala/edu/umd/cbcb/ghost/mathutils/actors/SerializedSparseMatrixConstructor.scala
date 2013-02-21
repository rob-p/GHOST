import org.apache.commons.math.linear._

import scala.actors.Actor

case class AddEntry( i: Int, j: Int, e: Double )
// Tells the actor that they are done
case class Stop()

class SerializedSparseSymmetricMatrixConstructor[MatT <: RealMatrix]( val m: MatT ) extends Actor { 
  

  def act() = { 
    
    loopWhile( !done ) { 
      react { 
	case AddEntry(i, j, e) => { 
	  m.setEntry(i,j,e)
	  m.setEntry(j,i,e)
	}
	case Stop() => { 
	  exit()
	}
      } // end react
    } // end loopWhile

  }
  
  def done() = { true }

}

