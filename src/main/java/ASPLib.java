package asp;

import com.sun.jna.Library;
import com.sun.jna.Native;
import com.sun.jna.Platform;
import com.sun.jna.Pointer;

// This is the standard, stable way of mapping, which supports extensive
// customization and mapping of Java to native types.
public interface ASPLib extends Library {
    ASPLib INSTANCE = (ASPLib)
	Native.loadLibrary("asp", ASPLib.class);
    
    void asp(int size, Pointer[] array, int[] col_mate, int[] row_mate);
}
