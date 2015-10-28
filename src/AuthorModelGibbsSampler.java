/*
 * AuthorModelGibbsSampler is free software; you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the Free
 * Software Foundation; either version 2 of the License, or (at your option) any
 * later version.
 *
 * AuthorModelGibbsSampler is distributed in the hope that it will be useful, but
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
 *
 *  Acknowledgement: Gregor Heinrich (gregor :: arbylon : net) 
 *  (This file is adapted from the Java code of Gregor Heinrich (gregor@arbylon.net)
 *    http://www.arbylon.net/projects/LdaGibbsSampler.java
 *  part of the org.knowceans experimental software packages.)
 *
 */
 
import java.text.DecimalFormat;
import java.text.NumberFormat;

import java.io.*;

/**
 *
 * Gibbs sampler for estimating the best assignments of topics for words, entities and
 * documents in a corpus. The model is discussed in Michal Rosen-Zvi' paper
 * "The Author-Topic Model for Authors and Documents"
 * (2004).
 *
 * @author mranahmd
 *
 */
public class AuthorModelGibbsSampler {

    /**
     * document data (term lists)
     */
    int[][] documentWords;
    
    /**
     * document data (entity lists)
     */
    int[][] documentAuthors;

    /**
     * vocabulary size of words
     */
    int W;
    
    /**
     * vocabulary size of Authors
     */
    int A;
    
    /**
     * Dirichlet parameter (author--word associations)
     */
    double beta;

    /**
     * author assignments for each word.
     */
    int x[][];
    
    /**
     * caw[i][j] number of instances of word i (term?) assigned to author j.
     */
    int[][] cwa;

    /**
     * sumcaw[j] total number of words assigned to author j.
     */
    int[] sumcwa;

    /**
     * cumulative statistics of phi
     */
    double[][] phisum;

    /**
     * size of statistics
     */
    int numstats;

    /**
     * sampling lag (?)
     */
    private static int THIN_INTERVAL = 20;

    /**
     * burn-in period
     */
    private static int BURN_IN = 100;

    /**
     * max iterations
     */
    private static int ITERATIONS = 1000;

    /**
     * sample lag (if -1 only one sample taken)
     */
    private static int SAMPLE_LAG;

    private static int dispcol = 0;

    /**
     * Initialise the Gibbs sampler with data.
     * 
     * @param W
     *            word vocabulary size
     * @param data
     *
	 * @param E
     *            entity vocabulary size
     * @param data
     */
    public AuthorModelGibbsSampler(int[][] documentWords, int[][] documentAuthors, int W, int A) {
        this.documentWords = documentWords;
        this.documentAuthors = documentAuthors;
        this.W = W;
        this.A = A;		
    }

    /**
     * Initialisation: Must start with an assignment of observations to topics ?
     * Many alternatives are possible, I chose to perform random assignments
     * with equal probabilities
     * 
     * @param K
     *            number of topics
     * @return z assignment of topics to words
     */
    public void initialState() {
        int i;

        int M = documentWords.length;

        // initialise count variables.
        cwa = new int[W][A];
        sumcwa = new int[A];
        
        // The x_i are are initialised to values in [1,A] to determine the
        // initial state of the Markov chain.

        x = new int[M][];
        for (int m = 0; m < M; m++) {
            int N = documentWords[m].length;
            x[m] = new int[N];
            for (int n = 0; n < N; n++) {
                int author = documentAuthors[m][(int) (Math.random() * documentAuthors[m].length)];
                x[m][n] = author;
                // number of instances of entity i assigned to topic j
                cwa[documentWords[m][n]][author]++;
                // total number of entities assigned to topic j.
                sumcwa[author]++;            
            }
        }
    }

    /**
     * Main method: Select initial state ? Repeat a large number of times: 1.
     * Select an element 2. Update conditional on other elements. If
     * appropriate, output summary for each run.
     * 
     * @param K
     *            number of topics
     * @param beta
     *            symmetric prior parameter on topic--word associations
     * @param alpha
     *            symmetric prior parameter on topic--author associations
     */
    private void gibbs(double beta) {
        this.beta = beta;

        // init sampler statistics
        if (SAMPLE_LAG > 0) {
            //phisum = new double[W][A];
            numstats = 0;
        }

        // initial state of the Markov chain:
        initialState();

        System.out.println("Sampling " + ITERATIONS
            + " iterations with burn-in of " + BURN_IN + " (B/S="
            + THIN_INTERVAL + ").");

        for (int i = 0; i < ITERATIONS; i++) {

            // for all z_i
            for (int m = 0; m < x.length; m++) {
                for (int n = 0; n < x[m].length; n++) {
                    int author = sampleFullConditional(m, n);
                    x[m][n] = author;	//word-author
                }
            }

            if ((i < BURN_IN) && (i % THIN_INTERVAL == 0)) {
                System.out.print("B");
                dispcol++;
            }
            // display progress
            if ((i > BURN_IN) && (i % THIN_INTERVAL == 0)) {
                System.out.print("S");
                dispcol++;
            }
            // get statistics after burn-in
            if ((i > BURN_IN) && (SAMPLE_LAG > 0) && (i % SAMPLE_LAG == 0)) {
                //updateParams();
                System.out.print("|");
                if (i % THIN_INTERVAL != 0)
                    dispcol++;
            }
            if (dispcol >= 100) {
                System.out.println();
                dispcol = 0;
            }
        }
    }
    
    /**
     * Sample an author x_i from the full conditional distribution: p(x_i = k |
     * x_-1, w) = (cwa_-i,j + beta)/(sumcwa- + W * beta)
     * 
     * @param m
     *            document
     * @param n
     *            word
     */
    private int sampleFullConditional(int m, int n) {
        
        // remove x_i from the count variables
        int author = x[m][n];
        cwa[documentWords[m][n]][author]--;
        sumcwa[author]--;
                
        // do multinomial sampling via cumulative method:
        double[] p = new double[documentAuthors[m].length];
        for (int k = 0; k < documentAuthors[m].length; k++) {
            p[k] = (cwa[documentWords[m][n]][documentAuthors[m][k]] + beta) / (sumcwa[documentAuthors[m][k]] + W * beta);
        }
        // cumulate multinomial parameters
        for (int k = 1; k < p.length; k++) {
            p[k] += p[k - 1];
        }
        // scaled sample because of unnormalised p[]
        double u = Math.random() * p[p.length - 1];
        for (int k = 0; k < p.length; k++) {
            if (u < p[k]){
            	author = documentAuthors[m][k];
                break;
            }
        }

        // add newly estimated x_i to count variables
        cwa[documentWords[m][n]][author]++;
        sumcwa[author]++;      

   	return author;
    }

    /**
     * Add to the statistics the values of theta and phi for the current state.
     */
    private void updateParams() {
        for (int w = 0; w < W; w++) {
        	for (int k = 0; k < A; k++) {
                phisum[w][k] += (cwa[w][k] + beta) / (sumcwa[k] + W * beta);
            }
        }
        numstats++;
    }

    /**
     * Retrieve estimated author--word associations. If sample lag > 0 then the
     * mean value of all sampled statistics for phi[][] is taken.
     * 
     * @return phi multinomial mixture of topic words (A x W)
     */
    public double[][] getPhi() {
        double[][] phi = new double[W][A];
        if (SAMPLE_LAG > 0) {
        	for (int w = 0; w < W; w++) {
        		for (int k = 0; k < A; k++) {
                    phi[w][k] = phisum[w][k] / numstats;
                }
            }
        } else {
       		for (int w = 0; w < W; w++) {
        		for (int k = 0; k < A; k++) {
                    phi[w][k] = (cwa[w][k] + beta) / (sumcwa[k] + W * beta);
                }
            }
        }
        return phi;
    }
    
     /**
     * 
     * @return cwa (K x A)
     */
    public int[][] getCwa() {
        return cwa;
    }

	/**
     * 
     * @return cwa (K x A)
     */
    public int[][] getCwaUrn() {
    	int[][] cwaUrn;
    	cwaUrn = new int[A][];
    	int[] lastIndex = new int[A];
    	for(int a=0; a<A;a++){
    		cwaUrn[a] = new int[sumcwa[a]];
    	}
    	int author, word;
        for (int m = 0; m < x.length; m++) {
            for (int n = 0; n < x[m].length; n++) {  
            	author = x[m][n];
            	word = documentWords[m][n];
            	cwaUrn[author][lastIndex[author]] = word;
            	lastIndex[author]++;
            }
        }  	
        return cwaUrn;
    }

    /**
     * Configure the gibbs sampler
     * 
     * @param iterations
     *            number of total iterations
     * @param burnIn
     *            number of burn-in iterations
     * @param thinInterval
     *            update statistics interval
     * @param sampleLag
     *            sample interval (-1 for just one sample at the end)
     */
    public void configure(int iterations, int burnIn, int thinInterval,
        int sampleLag) {
        ITERATIONS = iterations;
        BURN_IN = burnIn;
        THIN_INTERVAL = thinInterval;
        SAMPLE_LAG = sampleLag;
    }

    public static void main(String[] args) {
	if (args.length != 5)
        {
            System.err.println("\nUSAGE: java AuthorModelGibbsSampler <documentWordIndexFile> <wordVocabFile> <documentAuthorIndexFile> <authorVocabFile> <beta>\n");
            System.exit(0);
        }
        
        try
    	{
			CorpusReader wcr = new CorpusReader(args[0], args[1]);
			int W = wcr.getV();
			int M = wcr.getD();
			int[][] documentWords = wcr.getDocumentTermIndex();
		
			CorpusReader acr = new CorpusReader(args[2], args[3]);
			int A = acr.getV();
			int[][] documentAuthors = acr.getDocumentTermIndex();
			if(M != acr.getD()){
				System.out.println("\n Fatal Error: Mismatch between word documents size " + M + " and author document size " + acr.getD());
			}

			System.out.println("\n Fatal Error: UpdateParams Commented\n");

			double beta = Float.parseFloat(args[4]);

			System.out.println("Author Model using Gibbs Sampling.");
			System.out.println("\nTraining Author Models for " + M + " documents with " + W + " words vocabulary and " + A + " author vocabulary ...\n");

			AuthorModelGibbsSampler am = new AuthorModelGibbsSampler(documentWords, documentAuthors, W, A);
			am.configure(2500, 1000, 100, 10);
			am.gibbs(beta);

			double[][] phi = am.getPhi();
			//int[][] cwa = am.getCwa();

			//TopicUtils.saveCSV(phi,"./AM_phi");
			//TopicUtils.saveCSV(cwa,"./AM_cwa");
			TopicUtils.saveTopAuthorWords(phi, wcr.getVocab(), acr.getVocab(), "./AM_top-author-words", 25);

		}catch (IOException ex){
            System.err.println("\nIOException....\nUSAGE: java AuthorModelGibbsSampler <documentWordIndexFile> <wordVocabFile> <documentAuthorIndexFile> <authorVocabFile> <beta>\n");
			System.err.println(ex.getMessage());
		}
    }
}
