/*
 * AuthorTopicGibbsSampler is free software; you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the Free
 * Software Foundation; either version 2 of the License, or (at your option) any
 * later version.
 *
 * AuthorTopicGibbsSampler is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program; if not, write to the Free Software Foundation, Inc., 59 Temple
 * Place, Suite 330, Boston, MA 02111-1307 USA
 */

/*
 * Created on Aug 3, 2014
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
 * documents in a corpus. The algorithm is introduced in Michal Rosen-Zvi' paper
 * "The Author-Topic Model for Authors and Documents"
 * (2004).
 *
 * @author mranahmd
 * 
 */
public class AuthorTopicGibbsSampler {

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
     * number of topics
     */
    int K;
    

    /**
     * Dirichlet parameter (topic--word associations)
     */
    double beta;
    
    /**
     * Dirichlet parameter (topic--author associations)
     */
    double alpha;

    /**
     * topic assignments for each word.
     */
    int z[][];

    /**
     * author assignments for each word.
     */
    int x[][];

    /**
     * cwt[i][j] number of instances of word i (term?) assigned to topic j.
     */
    int[][] cwt;

    /**
     * sumcwt[j] total number of words assigned to topic j.
     */
    int[] sumcwt;
    
    /**
     * cat[i][j] number of instances of author i (term?) assigned to topic j.
     */
    int[][] cat;

    /**
     * sumcat[j] total number of author assigned to topic j.
     */
    int[] sumcat;

    /**
     * cumulative statistics of phi
     */
    double[][] phisum;
    
    /**
     * cumulative statistics of theta
     */
    double[][] thetasum;

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
    public AuthorTopicGibbsSampler(int[][] documentWords, int[][] documentAuthors, int W, int A) {
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
    public void initialState(int K) {
        int i;

        int M = documentWords.length;

        // initialise count variables.
        cwt = new int[W][K];
        cat = new int[A][K];
        sumcwt = new int[K];
        sumcat = new int[A];
        
        // The z_i are are initialised to values in [1,K] to determine the
        // initial state of the Markov chain.

        z = new int[M][];
        x = new int[M][];
        for (int m = 0; m < M; m++) {
            int N = documentWords[m].length;
            z[m] = new int[N];
            x[m] = new int[N];
            for (int n = 0; n < N; n++) {
                int topic = (int) (Math.random() * K);
                z[m][n] = topic;
                // number of instances of word i assigned to topic j
                cwt[documentWords[m][n]][topic]++;
                // total number of words assigned to topic j.
                sumcwt[topic]++;
                
                int author = documentAuthors[m][(int) (Math.random() * documentAuthors[m].length)];
                x[m][n] = author;
                // number of instances of entity i assigned to topic j
                cat[author][topic]++;
                // total number of entities assigned to topic j.
                sumcat[author]++;            
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
    private void gibbs(int K, double alpha, double beta) {
        this.K = K;
        this.beta = beta;
        this.alpha = alpha;

        // init sampler statistics
        if (SAMPLE_LAG > 0) {
            phisum = new double[K][W];
            thetasum = new double[K][A];
            numstats = 0;
        }

        // initial state of the Markov chain:
        initialState(K);

        System.out.println("Sampling " + ITERATIONS
            + " iterations with burn-in of " + BURN_IN + " (B/S="
            + THIN_INTERVAL + ").");

        for (int i = 0; i < ITERATIONS; i++) {

            // for all z_i
            for (int m = 0; m < z.length; m++) {
                for (int n = 0; n < z[m].length; n++) {
                    int at[] = sampleFullConditional(m, n);
                    z[m][n] = at[0];	//word-topic
                    x[m][n] = at[1];	//word-author
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
                updateParams();
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
     * Sample a topic z_i and x_i from the full conditional distribution: p(z_i = j, x_i = k |
     * x_-1, w, z_-i) = (cat_-i,k + alpha)/(sumcat_-i,j(.) + K * alpha) * (cwt_-i,j + beta)/(sumcwt- + W * beta)
     * 
     * @param m
     *            document
     * @param n
     *            word
     */
    private int[] sampleFullConditional(int m, int n) {

        // remove z_i from the count variables
        int topic = z[m][n];
        cwt[documentWords[m][n]][topic]--;
        sumcwt[topic]--;
        
        // remove x_i from the count variables
        int author = x[m][n];
        cat[author][topic]--;
        sumcat[author]--;
                
        // do multinomial sampling via cumulative method:
        double[][] p = new double[K][documentAuthors[m].length];
        for (int j = 0; j < K; j++) {
        	for (int k = 0; k < documentAuthors[m].length; k++) {
            	p[j][k] = (cat[documentAuthors[m][k]][j] + alpha) / (sumcat[documentAuthors[m][k]] + K * alpha)
                		  * (cwt[documentWords[m][n]][j] + beta) / (sumcwt[j] + W * beta);
            }
        }
        // cumulate multinomial parameters
        double cump = 0;
        for (int j = 0; j < K; j++) {
        	for (int k = 0; k < documentAuthors[m].length; k++) {            
        		p[j][k] += cump;
        		cump = p[j][k];
        	}
        }
        // scaled sample because of unnormalised p[]
        double u = Math.random() * p[K - 1][documentAuthors[m].length - 1];
        boolean breakFlag = false;
        for (int j = 0; j < K; j++) {
        	for (int k = 0; k < documentAuthors[m].length; k++) {            
        		if (u < p[j][k]){
        			author = documentAuthors[m][k];
        			topic = j;
		            breakFlag = true;
		            break;
		        }
            }
            if (breakFlag){
	            break;
		    }
        }

        // add newly estimated ze_i to count variables
        cwt[documentWords[m][n]][topic]++;
        sumcwt[topic]++;
        cat[author][topic]++;
        sumcat[author]++;        

   		return new int[] {topic, author};
    }

    /**
     * Add to the statistics the values of theta and phi for the current state.
     */
    private void updateParams() {
        for (int k = 0; k < K; k++) {
            for (int a = 0; a < A; a++) {
                thetasum[k][a] += (cat[a][k] + alpha) / (sumcat[a] + K * alpha);
            }
        }
        for (int k = 0; k < K; k++) {
            for (int w = 0; w < W; w++) {
                phisum[k][w] += (cwt[w][k] + beta) / (sumcwt[k] + W * beta);
            }
        }
        numstats++;
    }

    /**
     * Retrieve estimated document--topic associations. If sample lag > 0 then
     * the mean value of all sampled statistics for theta[][] is taken.
     * 
     * @return theta multinomial mixture of document topics (A x K)
     */
    public double[][] getTheta() {
        double[][] theta = new double[K][A];

        if (SAMPLE_LAG > 0) {
			for (int k = 0; k < K; k++) {
            	for (int a = 0; a < A; a++) {
                    theta[k][a] = thetasum[k][a] / numstats;
                }
            }

        } else {
			for (int k = 0; k < K; k++) {
            	for (int a = 0; a < A; a++) {
                    theta[k][a] = (cat[a][k] + alpha) / (sumcat[a] + K * alpha);
                }
            }
        }

        return theta;
    }

    /**
     * Retrieve estimated topic--word associations. If sample lag > 0 then the
     * mean value of all sampled statistics for phi[][] is taken.
     * 
     * @return phi multinomial mixture of topic words (K x W)
     */
    public double[][] getPhi() {
        double[][] phi = new double[K][W];
        if (SAMPLE_LAG > 0) {
            for (int k = 0; k < K; k++) {
                for (int w = 0; w < W; w++) {
                    phi[k][w] = phisum[k][w] / numstats;
                }
            }
        } else {
            for (int k = 0; k < K; k++) {
                for (int w = 0; w < W; w++) {
                    phi[k][w] = (cwt[w][k] + beta) / (sumcwt[k] + W * beta);
                }
            }
        }
        return phi;
    }
    
     /**
     * 
     * @return cwt (V x K)
     */
    public int[][] getCwt() {
        return cwt;
    }
    
     /**
     * 
     * @return cet (A x K)
     */
    public int[][] getCat() {
        return cat;
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
	if (args.length != 7)
        {
            System.err.println("\nUSAGE: java ATGibbsSampler <documentWordIndexFile> <wordVocabFile> <documentAuthorIndexFile> <authorVocabFile> <K> <alpha> <beta>\n");
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

			int K = Integer.parseInt(args[4]);
			double alpha = Float.parseFloat(args[5]);
			double beta = Float.parseFloat(args[6]);

			System.out.println("Author-Topic Model using Gibbs Sampling.");
			System.out.println("\nTraining " + K + " AT topics for " + M + " documents with " + W + " words vocabulary and " + A + " author vocabulary ...\n");

			AuthorTopicGibbsSampler atlda = new AuthorTopicGibbsSampler(documentWords, documentAuthors, W, A);
			atlda.configure(2500, 1000, 100, 10);
			atlda.gibbs(K, alpha, beta);

			double[][] theta = atlda.getTheta();
			double[][] phi = atlda.getPhi();
			int[][] cwt = atlda.getCwt();
			int[][] cat = atlda.getCat();

			TopicUtils.saveCSV(theta,"./AT_theta");
			TopicUtils.saveCSV(phi,"./AT_phi");
			TopicUtils.saveCSV(cwt,"./AT_cwt");
			TopicUtils.saveCSV(cat,"./AT_cat");
			TopicUtils.saveTopTopicWords(phi, wcr.getVocab(), "./AT_top-topic-words", 10);
			TopicUtils.saveTopTopicWords(theta, acr.getVocab(), "./AT_top-topic-authors", 10);

		}catch (IOException ex){
            System.err.println("\nIOException....\nUSAGE: java ATGibbsSampler <documentWordIndexFile> <wordVocabFile> <documentAuthorIndexFile> <authorVocabFile> <K> <alpha> <beta>\n");
			System.err.println(ex.getMessage());
		}
    }
}
