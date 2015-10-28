/*
 * ArrayMatrixUtils.java is free software; you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the Free
 * Software Foundation; either version 2 of the License, or (at your option) any
 * later version.
 *
 * ArrayMatrixUtils.java is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program; if not, write to the Free Software Foundation, Inc., 59 Temple
 * Place, Suite 330, Boston, MA 02111-1307 USA
 */

/*
 *  Created on Aug 3, 2014
 *
 * 	Imran SHEIKH
 *	imran.sheikh@loria.fr, mranahmd@gmail.com
 * 	Multispeech, LORIA-INRIA
 * 	Nancy, France
 */

import java.util.*;
import java.io.*;

/**
 *
 * Utility functions for handling arrays and matrices for topic models in this repository.
 *
 * @author mranahmd
 *
 */
public class ArrayMatrixUtils {
	public static double[][] transposeMatrix(double [][] m){
        double[][] temp = new double[m[0].length][m.length];
        for (int i = 0; i < m.length; i++)
            for (int j = 0; j < m[0].length; j++)
                temp[j][i] = m[i][j];
        return temp;
    }

	public static double[][] matadd(double [][] m, double [][] n){
        double[][] temp = new double[m.length][m[0].length];
	if((m.length != n.length) || (m[0].length != n[0].length)){
		System.err.println("\n Matrix Addition : Fatal Error: Mismatch between matrix dimensions [" + m.length + "x" + m[0].length + "] [" + n.length + "x" + n[0].length + "]\n");
		return temp;
	}
        for (int i = 0; i < m.length; i++)
            for (int j = 0; j < m[0].length; j++)
                temp[i][j] = m[i][j] + n[i][j];
        return temp;
    }

	public static double[] getVectorMagnitude(double[][] mat){
		double[] ret = new double[mat.length];
		double sumsqr = 0;
    	for(int i=0; i< mat.length; i++){
    		sumsqr = 0;
    		for(int j=0; j< mat[i].length; j++){
    			sumsqr += mat[i][j] * mat[i][j];
    		}	
    		ret[i] = Math.sqrt(sumsqr);
		}
		return ret;
	}

	public static int[][] transposeMatrix(int [][] m){
        int[][] temp = new int[m[0].length][m.length];
        for (int i = 0; i < m.length; i++)
            for (int j = 0; j < m[0].length; j++)
                temp[j][i] = m[i][j];
        return temp;
    }

    public static double sumOfMatrixElements(double [][] m){
        double sum = 0;
        for (int i = 0; i < m.length; i++)
            for (int j = 0; j < m[0].length; j++)
                sum += m[i][j];
        return sum;
    }
	
	public static void SortWithIndex(double[] arr, Integer[] i){
		final Integer[] idx = i;
		final double[] data = arr;

		Arrays.sort(idx, new Comparator<Integer>() {
    		@Override public int compare(final Integer o1, final Integer o2) {
        	return Double.compare(data[o1], data[o2]);
    	}
		});
	}
	
	public static double[][] readMatrixFromFile(String filename) throws IOException{

		int numRows = countLines(filename);
		double[][] mat = new double[numRows][];
		int r = 0;
		BufferedReader f = new BufferedReader(new FileReader(filename));
     		String line = f.readLine();
			while(line != null){
				line = line.replaceAll("(\\r|\\n)", "");	
				line = line.trim();	
				if (line.isEmpty() || line.trim().equals("")){
					System.err.println("Fatal Error: Blank line in " + filename);			
				}
				String [] cols = line.split("\\s+");
				mat[r] = new double[cols.length];
				for (int c = 0; c < cols.length; c++) {
						mat[r][c] = Double.parseDouble(cols[c]);
				}
			r++; 
			line = f.readLine();
		}
		return mat;
	}
		
	public static int[][] readIntMatrixFromFile(String filename) throws IOException{

		int numRows = countLines(filename);
		int[][] mat = new int[numRows][];
		int r = 0;
		BufferedReader f = new BufferedReader(new FileReader(filename));
     		String line = f.readLine();
			while(line != null){
				line = line.replaceAll("(\\r|\\n)", "");	
				line = line.trim();	
				if (line.isEmpty() || line.trim().equals("")){
					System.err.println("Fatal Error: Blank line in " + filename);			
				}
				String [] cols = line.split("\\s+");
				mat[r] = new int[cols.length];
				for (int c = 0; c < cols.length; c++) {
						mat[r][c] = Integer.parseInt(cols[c]);
				}
			r++; 
			line = f.readLine();
		}
		return mat;
	}	
	
	public static int countLines(String filename) throws IOException {
		BufferedReader br = new BufferedReader(new FileReader(filename));
		String line;
		int lines = 0;
		while ((line = br.readLine()) != null) {
			lines++;
  			if (line.trim().isEmpty() || line.trim().equals("")){
   				System.err.println("Fatal Error: Blank line in " + filename + " at " + lines);			
  			}
		}
		return lines;
	}
	
	public int countLinesFast(String filename) throws IOException {
   		InputStream is = new BufferedInputStream(new FileInputStream(filename));
    	try {
        	byte[] c = new byte[1024];
        	int count = 0;
        	int readChars = 0;
        	boolean endsWithoutNewLine = false;
        	while ((readChars = is.read(c)) != -1) {
            	for (int i = 0; i < readChars; ++i) {
                	if (c[i] == '\n')
                    	++count;
            	}
            	endsWithoutNewLine = (c[readChars - 1] != '\n');
        	}
        	if(endsWithoutNewLine) {
            	++count;
        	} 
        	return count;
    	} finally {
        	is.close();
    	}
	}
	
	public static double[][] matmul(double[][] a, double[][] b) {
		if(a[0].length != b.length){
			System.err.println("\n Matrix Multiplication : Fatal Error: Mismatch between matrix dimensions [" + a.length + "x" + a[0].length + "] [" + b.length + "x" + b[0].length + "]\n");
		}
		int m = a.length, n = a[0].length, p = b[0].length;
		double[][] x = new double[m][p];
		double[][] c = new double[p][n];
		for (int i = 0; i < n; ++i) // transpose
			for (int j = 0; j < p; ++j)
				c[j][i] = b[i][j];
		for (int i = 0; i < m; ++i)
			for (int j = 0; j < p; ++j) {
				double s = 0.0;
				for (int k = 0; k < n; ++k)
					s += a[i][k] * c[j][k];
				x[i][j] = s;
			}
		return x;
	}
	
	public static int[][] matmul(int[][] a, int[][] b) {
		if(a[0].length != b.length){
			System.out.println("\n Fatal Error: Mismatch between matrix dimensions [" + a.length + "x" + a[0].length + "] [" + b.length + "x" + b[0].length + "]\n");
		}
		int m = a.length, n = a[0].length, p = b[0].length;
		int[][] x = new int[m][p];
		int[][] c = new int[p][n];
		for (int i = 0; i < n; ++i) // transpose
			for (int j = 0; j < p; ++j)
				c[j][i] = b[i][j];
		for (int i = 0; i < m; ++i)
			for (int j = 0; j < p; ++j) {
				int s = 0;
				for (int k = 0; k < n; ++k)
					s += a[i][k] * c[j][k];
				x[i][j] = s;
			}
		return x;
	}
	
	
	/**
     * 
     * @return cwaUrn (K x A)
     */
    public static void saveCwaUrn(int[][] cwa, String filename)
    	 throws IOException
    {
        PrintStream fileStream = new PrintStream(filename);
        for( int a = 0; a < cwa.length; a++ ) {
            for( int w = 0; w < cwa[a].length; w++ ) {
            	while(cwa[a][w] != 0){
                	fileStream.print(w+" ");
                	cwa[a][w]--;
                }
            }
            fileStream.println();
        }
        fileStream.close();   
    }
	
	/**
     * Saves a double matrix to disk using in a Column Space Value (CSV) format.
     *
     * @param A The matrix being saved.
     * @param fileName Name of the file its being saved at.
     * @throws java.io.IOException
     */
    public static void saveCSV(double[][] A , String fileName )
        throws IOException
    {
    	int cols = A[0].length;
    	boolean isMatrix = true;
        PrintStream fileStream = new PrintStream(fileName);
        for( int i = 0; i < A.length; i++ ) {
        	if(A[i].length != cols){isMatrix = false;}
            for( int j = 0; j < A[i].length; j++ ) {
                fileStream.print(A[i][j]+" ");
            }
            fileStream.println();
        }
        fileStream.close();
        if(isMatrix){
        	System.out.println("\nSaved " + A.length + "x" + A[0].length + " matrix to "+ fileName);
        }else{
        	System.out.println("\nSaved " + A.length + " row 2d array to "+ fileName);
        }
    }
    
    /**
     * Saves an int matrix to disk using in a Column Space Value (CSV) format. 
     *
     * @param A The matrix being saved.
     * @param fileName Name of the file its being saved at.
     * @throws java.io.IOException
     */
    public static void saveCSV(int[][] A , String fileName )
        throws IOException
    {
    	int cols = A[0].length;
    	boolean isMatrix = true;
        PrintStream fileStream = new PrintStream(fileName);
        for( int i = 0; i < A.length; i++ ) {
        	if(A[i].length != cols){isMatrix = false;}
            for( int j = 0; j < A[i].length; j++ ) {
                fileStream.print(A[i][j]+" ");
            }
            fileStream.println();
        }
        fileStream.close();
        if(isMatrix){
        	System.out.println("\nSaved " + A.length + "x" + A[0].length + " matrix to "+ fileName);
        }else{
        	System.out.println("\nSaved " + A.length + " row 2d array to "+ fileName);
        }
    }

  /**
  * A dummy main to test matrix multiplication	
  **/
  public static void main(String[] args) {
	if (args.length != 6)
        {
            System.err.println("\nUSAGE: java ArrayMatrixUtils <matFile1> <matFile2> <matFile1-TransposeFlag_0|1> <matFile2-TransposeFlag_0|1> <double-0|Int-1_Flag> <op-tag>\n");
            System.exit(0);
        }
        
        try
    	{
			
			if(Integer.parseInt(args[4]) == 1){
				int[][] A, B;
				if(Integer.parseInt(args[2]) == 0){
					A = ArrayMatrixUtils.readIntMatrixFromFile(args[0]);
				}else{
					A = ArrayMatrixUtils.transposeMatrix(ArrayMatrixUtils.readIntMatrixFromFile(args[0]));
				}
				if(Integer.parseInt(args[3]) == 0){
					B = ArrayMatrixUtils.readIntMatrixFromFile(args[1]);
				}else{
					B = ArrayMatrixUtils.transposeMatrix(ArrayMatrixUtils.readIntMatrixFromFile(args[1]));
				}
				int[][] cwa = ArrayMatrixUtils.matmul(A, B);
				saveCSV(cwa, args[5]);
			}else{
				double[][] A, B;
				if(Integer.parseInt(args[2]) == 0){
					A = ArrayMatrixUtils.readMatrixFromFile(args[0]);
				}else{
					A = ArrayMatrixUtils.transposeMatrix(ArrayMatrixUtils.readMatrixFromFile(args[0]));
				}
				if(Integer.parseInt(args[3]) == 0){
					B = ArrayMatrixUtils.readMatrixFromFile(args[1]);
				}else{
					B = ArrayMatrixUtils.transposeMatrix(ArrayMatrixUtils.readMatrixFromFile(args[1]));
				}
				saveCSV(ArrayMatrixUtils.matmul(A, B), args[5]);
			}
		}catch (IOException ex){
            System.err.println("\nUSAGE: java ArrayMatrixUtils <matFile1> <matFile2> <matFile1-TransposeFlag_0|1> <matFile2-TransposeFlag_0|1> <double-0|Int-1_Flag> <op-tag>\n");
            System.err.println(ex.getMessage());
		}
    }
}
