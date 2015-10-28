/*
 * TopicUtils.java is free software; you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the Free
 * Software Foundation; either version 2 of the License, or (at your option) any
 * later version.
 *
 * TopicUtils.java is distributed in the hope that it will be useful, but
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

import java.io.*;
import java.util.Arrays;

/**
 * Utility functions for learning and checking topic models in this repository.
 * 
 * @author mranahmd
 *
 */
public class TopicUtils {    
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
     * Saves a double vector to disk using in a Column Space Value (CSV) format. 
     *
     * @param A The matrix being saved.
     * @param fileName Name of the file its being saved at.
     * @throws java.io.IOException
     */
    public static void saveCSV(double[] A , String fileName )
        throws IOException
    {
        PrintStream fileStream = new PrintStream(fileName);
        for( int i = 0; i < A.length; i++ ) {
            fileStream.print(A[i] + " ");
        }
        fileStream.close();
        System.out.println("\nSaved " + A.length + " row vector to "+ fileName);
    }
    
    /**
     * Saves an int vector to disk using in a Column Space Value (CSV) format. 
     *
     * @param A The matrix being saved.
     * @param fileName Name of the file its being saved at.
     * @throws java.io.IOException
     */
    public static void saveCSV(int[] A , String fileName )
        throws IOException
    {
        PrintStream fileStream = new PrintStream(fileName);
        for( int i = 0; i < A.length; i++ ) {
            fileStream.print(A[i] + " ");
        }
        fileStream.close();
        System.out.println("\nSaved " + A.length + " row vector to "+ fileName);
    }
    
    /**
     * Saves top words for each topic. Helps to check the topics learnt. 
     *
     * @param phi The phi matrix.
     * @param vocab The word/entity vocabulary file
     * @param fileName Name of the file its being saved at.
     * @param numWordsPerTopic Number of top words to save per topic 
     * @throws java.io.IOException
     */    
     public static void saveTopTopicWords(double[][] phi, String[] vocab, String fileName, int numWordsPerTopic)
        throws IOException
    {
    	double[][] A = phi;
    
        Integer[] idx = new Integer[A[0].length];
      	for( int i = 0; i < idx.length; i++ ) {
			idx[i] = i;
		}
		
        PrintStream fileStream = new PrintStream(fileName);
        for( int i = 0; i < A.length; i++ ) {
			Integer[] idxtmp = idx;
			ArrayMatrixUtils.SortWithIndex(A[i], idxtmp);
            for( int j = 0; j< numWordsPerTopic; j++ ) {
                fileStream.print(vocab[idxtmp[idxtmp.length - 1 - j]]+" ");
            }
            fileStream.println();
        }
        fileStream.close();
        
    	System.out.println("\nSaved top " + numWordsPerTopic + " terms for each of the "+ phi.length + " topics to "+ fileName);
    }
    
    /**
     * Saves top topics each supertopic in CorrLda2. Helps to check the supertopics learnt. 
     *
     * @param psi The psi matrix in CorrLda2.
     * @param fileName Name of the file its being saved at.
     * @param numTopicsPerSupertopic Number of top topics to save per supertopic 
     * @throws java.io.IOException
     */     public static void saveTopTopicsOfSupertopics(double[][] psi, String fileName, int numTopicsPerSupertopic)
        throws IOException
    {
    	double[][] A = psi;
    
        Integer[] idx = new Integer[A[0].length];
      	for( int i = 0; i < idx.length; i++ ) {
			idx[i] = i;
		}
		
        PrintStream fileStream = new PrintStream(fileName);
        for( int i = 0; i < A.length; i++ ) {
			Integer[] idxSorted = idx;
			ArrayMatrixUtils.SortWithIndex(A[i], idxSorted);
            for( int j = 0; j< numTopicsPerSupertopic; j++ ) {
                fileStream.print("T" +idxSorted[idx.length - 1 - j]+" ");
            }
            fileStream.println();
        }
        fileStream.close();
        System.out.println("\nSaved top " + numTopicsPerSupertopic + " topics for each of the "+ psi.length + " supertopics to "+ fileName);
    }    
    
    /**
     * Saves top words for each author. Helps to check the author model. 
     *
     * @param phi The phi matrix for author model.
     * @param wvocab The word vocabulary file
     * @param avocab The author names vocabulary file
     * @param fileName Name of the file its being saved at.
     * @param numWordsPerAuthor Number of top words to save per author 
     * @throws java.io.IOException
     */  
	public static void saveTopAuthorWords(double[][] phi, String[] wVocab, String[] aVocab, String fileName, int numWordsPerAuthor)
        throws IOException
    {
    	double[][] A = ArrayMatrixUtils.transposeMatrix(phi);	// -> (AxW), phi -> (WxA)
    
        Integer[] idx = new Integer[A[0].length];
      	for( int i = 0; i < idx.length; i++ ) {
			idx[i] = i;
		}
		
        PrintStream fileStream = new PrintStream(fileName);
        for( int i = 0; i < A.length; i++ ) {
			Integer[] idxtmp = idx;
			ArrayMatrixUtils.SortWithIndex(A[i], idxtmp);
			fileStream.print("[" + aVocab[i] + "] : ");
            for( int j = 0; j< numWordsPerAuthor; j++ ) {
                fileStream.print(wVocab[idxtmp[idxtmp.length - 1 - j]]+" ");
            }
            fileStream.println();
        }
        fileStream.close();
        
    	System.out.println("\nSaved top " + numWordsPerAuthor + " terms for each of the "+ A.length + " authors to "+ fileName);
    }
    
     /**
     * Calculate similarities between two set of document topic vectors. 
     *
     * @param thetaSetA One set of document topic vectors.
     * @param thetaSetB Second set of document topic vectors.
     */     
    public static double[][] calcDocumentCosineSimilarities(double[][] thetaSetA, double[][] thetaSetB){
    	double[][] ret;
    	if(thetaSetA[0].length != thetaSetB[0].length){
    		System.err.println("Fatal Error: Topic Dimensions do not agree ... ");
    		ret = new double[1][1];
    	}else{
    		int nA = thetaSetA.length;
    		int nB = thetaSetB.length;
    		double[][] thetaSetBT = ArrayMatrixUtils.transposeMatrix(thetaSetB);
    		ret = ArrayMatrixUtils.matmul(thetaSetA, thetaSetBT);
    		double setAMag [] = ArrayMatrixUtils.getVectorMagnitude(thetaSetA);
    		double setBMag [] = ArrayMatrixUtils.getVectorMagnitude(thetaSetB);
    		for(int i=0; i< nA; i++){
    			for(int j=0; j< nB; j++){
    				ret[i][j] = ret[i][j]/(setAMag[i] * setBMag[j]);
    			}
    		}
    	}
    	return ret;
    }
    
     /**
     * Read a vocabulary file into a String array. 
     *
     * @param filename Name of the vocabulary file.
     * @throws java.io.IOException
     */ 
    public static String[]readVocab(String filename) throws IOException{
		// Read the vocab file
		int V = countLines(filename);
		String[] dict = new String[V];
     	int v=0;
		BufferedReader vocab = new BufferedReader(new FileReader(filename));
     	String line = vocab.readLine();
		while(line != null){
		        line = line.replaceAll("(\\r|\\n)", "");	
		        line = line.trim();	
				if (line.isEmpty() || line.trim().equals("")){
					System.err.println("Fatal Error: Blank line in " + filename);
				}
		        dict[v] = line;
		        v++;
	       		line = vocab.readLine();
		}
		return dict;
	}
	
	/**
     * Count the lines in a file. 
     *
     * @param filename Name of the file.
     * @throws java.io.IOException
     */ 
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
}
