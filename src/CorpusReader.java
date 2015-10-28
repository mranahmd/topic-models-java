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

/**
 * Corpus related utility functions for learning topic models in this repository.
 * 
 * @author mranahmd
 *
 */
public class CorpusReader {
	int[][] documentTermIndex;
	int V, D;
	String[] dict;
	
    /**
     * Reads a corpus text file in document term-index format.
     * Each line in the corpus text file is a document with integer indexes for terms in the document.
     * Term indexes are zero-based and as per the specified vocabulary file (one word per line).
     *
     * @param documentTermIndexFile The corpus text file in term-index format.
     * @param vocabFile Name of the vocabulary file.
     * @throws java.io.IOException
     */	
     public CorpusReader(String documentTermIndexFile, String vocabFile) throws IOException{
		// Read the vocab file
		V = countLines(vocabFile);
		dict = new String[V];
     	int v=0;
		BufferedReader vocab = new BufferedReader(new FileReader(vocabFile));
     	String line = vocab.readLine();
		while(line != null){
		        line = line.replaceAll("(\\r|\\n)", "");	
		        line = line.trim();	
				if (line.isEmpty() || line.trim().equals("")){
					System.err.println("Fatal Error: Blank line in " + vocabFile);
				}
		        dict[v] = line;
		        v++;
	       		line = vocab.readLine();
		}

		// Read the documentTermIndexFile file
		D = countLines(documentTermIndexFile);
		documentTermIndex = new int[D][];
		int docNum = 0;
		BufferedReader terms = new BufferedReader(new FileReader(documentTermIndexFile));
     	line = terms.readLine();
		while(line != null){
			line = line.replaceAll("(\\r|\\n)", "");
			line = line.trim();	
			if (line.isEmpty() || line.trim().equals("")){
				System.err.println("Fatal Error: Blank line in " + documentTermIndexFile);
			}
			String [] lineWords = line.split("\\s+");
			documentTermIndex[docNum] = new int[lineWords.length];
			for (int i = 0; i < lineWords.length; i++) {
				documentTermIndex[docNum][i] = Integer.parseInt(lineWords[i]);
			}
			docNum++; 
			line = terms.readLine();
		}
	}
	
	/**
     * Count the lines in a file. 
     *
     * @param filename Name of the file.
     * @throws java.io.IOException
     */ 
	public int countLines(String filename) throws IOException {
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
	
	/**
     * Count the lines in a file. (Faster version). 
     *
     * @param filename Name of the file.
     * @throws java.io.IOException
     */ 
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
	
	public int[][] getDocumentTermIndex(){
		return documentTermIndex;
	}
	
	public int getV(){
		return V;
	}
	
	public int getD(){
		return D;
	}
	
	public String[] getVocab(){
		return dict;
	}
}


