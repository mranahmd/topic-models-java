/*
 * LdaGibbsSamplingForInference is free software; you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the Free
 * Software Foundation; either version 2 of the License, or (at your option) any
 * later version.
 *
 * LdaGibbsSamplingForInference is distributed in the hope that it will be useful, but
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
 *	imran.sheikh@loria.fr or mranahmd@gmail.com
 * 	Parole - LORIA
 * 	Nancy, France
 *
 *  Gregor Heinrich (gregor :: arbylon : net) 
 *  (This file is adapted from the Java code of Gregor Heinrich (gregor@arbylon.net)
 *    http://www.arbylon.net/projects/LdaGibbsSampler.java
 *  part of the org.knowceans experimental software packages.)
 */

import java.io.*;

/**
 * Gibbs sampler for infering the best assignments of topics for words and
 * documents in a test corpus. The algorithm is introduced in Tom Griffiths' paper
 * "Gibbs sampling in the generative model of Latent Dirichlet Allocation"
 * (2002).
 * 
 * @author mranahmd
 */
public class LdaGibbsSamplingForInference {

    /**
     * document data (term lists)
     */
    int[][] documents;

    /**
     * vocabulary size
     */
    int V;

    /**
     * number of topics
     */
    int K;

    /**
     * Dirichlet parameter (document--topic associations)
     */
    double alpha;

    /**
     * Dirichlet parameter (topic--term associations)
     */
    double beta;

    /**
     * topic assignments for each word.
     */
    int z[][];

    /**
     * nw[m][i][j] number of instances of word i (term?) assigned to topic j.
     */
    int[][][] nw;
    
    int[][] nwk;

    /**
     * cwt[i][j] number of instances of word i (term?) assigned to topic j in trained model.
     */
    int[][] cwt;
    
    /**
     * na[i][j] number of words in document i assigned to topic j.
     */
    int[][] nd;

    /**
     * nwsum[j] total number of words assigned to topic j in trained model.

     */
    int[] nwsum;

    /**
     * cwtsum[j] total number of words assigned to topic j in trained model.

     */
    int[] cwtsum;

    /**
     * nasum[i] total number of words in document i.
     */
    int[] ndsum;

    /**
     * cumulative statistics of theta
     */
    double[][] thetasum;

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
     * @param V
     *            vocabulary size
     * @param data
     */
    public LdaGibbsSamplingForInference(int[][] documents, int V) {

        this.documents = documents;
        this.V = V;
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

        int M = documents.length;

        // initialise count variables.
        nw = new int[M][V][K];
        nwk = new int[V][K];
        nd = new int[M][K];
        ndsum = new int[M];
        cwtsum = new int[K];
        nwsum = new int[K];

        // The z_i are are initialised to values in [1,K] to determine the
        // initial state of the Markov chain.

        z = new int[M][];
        for (int m = 0; m < M; m++) {
            int N = documents[m].length;
            z[m] = new int[N];
            for (int n = 0; n < N; n++) {
                int topic = (int) (Math.random() * K);
                z[m][n] = topic;
                // number of instances of word i assigned to topic j
                nw[m][documents[m][n]][topic]++;
                nwk[documents[m][n]][topic]++;
                // number of words in document i assigned to topic j.
                nd[m][topic]++;
                nwsum[topic]++;
            }
            // total number of words in document i
            ndsum[m] = N;
        }
        
        for (int k = 0; k < K; k++) {
            for (int w = 0; w < V; w++) {
                cwtsum[k] += cwt[w][k];
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
     * @param alpha
     *            symmetric prior parameter on document--topic associations
     * @param beta
     *            symmetric prior parameter on topic--term associations
     */
    private void gibbs(int K, double alpha, double beta, int[][] cwt) {
        this.K = K;
        this.alpha = alpha;
        this.beta = beta;
        this.cwt = cwt;
        
        // init sampler statistics
        if (SAMPLE_LAG > 0) {
            thetasum = new double[documents.length][K];
            phisum = new double[K][V];
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

                    // (z_i = z[m][n])
                    // sample from p(z_i|z_-i, w)
                    int topic = sampleFullConditional(m, n);
                    z[m][n] = topic;
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
     * Sample a topic z_i from the full conditional distribution: p(z_i = j |
     * z_-i, w) = (n_-i,j(w_i) + beta)/(n_-i,j(.) + W * beta) * (n_-i,j(d_i) +
     * alpha)/(n_-i,.(d_i) + K * alpha)
     * 
     * @param m
     *            document
     * @param n
     *            word
     */
    private int sampleFullConditional(int m, int n) {

        // remove z_i from the count variables
        // cwt is unchanged
        int topic = z[m][n];
        nw[m][documents[m][n]][topic]--;
        nwk[documents[m][n]][topic]--;
        nd[m][topic]--;
        ndsum[m]--;
	nwsum[topic]--;
	
        // do multinomial sampling via cumulative method:
        double[] p = new double[K];
        for (int k = 0; k < K; k++) {
            p[k] = (cwt[documents[m][n]][k] + nw[m][documents[m][n]][k] + beta) / (cwtsum[k] + nd[m][k] + V * beta)
                * (nd[m][k] + alpha) / (ndsum[m] + K * alpha);
        }
        // cumulate multinomial parameters
        for (int k = 1; k < p.length; k++) {
            p[k] += p[k - 1];
        }
        // scaled sample because of unnormalised p[]
        double u = Math.random() * p[K - 1];
        for (topic = 0; topic < p.length; topic++) {
            if (u < p[topic])
                break;
        }

        // add newly estimated z_i to count variables
        nw[m][documents[m][n]][topic]++;
        nwk[documents[m][n]][topic]++;
        nd[m][topic]++;
        ndsum[m]++;
	nwsum[topic]++;
        return topic;
    }

    /**
     * Add to the statistics the values of theta and phi for the current state.
     */
    private void updateParams() {
        for (int m = 0; m < documents.length; m++) {
            for (int k = 0; k < K; k++) {
                thetasum[m][k] += (nd[m][k] + alpha) / (ndsum[m] + K * alpha);
            }
        }
        for (int k = 0; k < K; k++) {
            for (int w = 0; w < V; w++) {
                phisum[k][w] += (cwt[w][k] + nwk[w][k] + beta) / (cwtsum[k] + nwsum[k] + V * beta);
            }
        }
        numstats++;
    }

    /**
     * Retrieve estimated document--topic associations. If sample lag > 0 then
     * the mean value of all sampled statistics for theta[][] is taken.
     * 
     * @return theta multinomial mixture of document topics (M x K)
     */
    public double[][] getTheta() {
        double[][] theta = new double[documents.length][K];

        if (SAMPLE_LAG > 0) {
            for (int m = 0; m < documents.length; m++) {
                for (int k = 0; k < K; k++) {
                    theta[m][k] = thetasum[m][k] / numstats;
                }
            }

        } else {
            for (int m = 0; m < documents.length; m++) {
                for (int k = 0; k < K; k++) {
                    theta[m][k] = (nd[m][k] + alpha) / (ndsum[m] + K * alpha);
                }
            }
        }

        return theta;
    }

    /**
     * Retrieve estimated topic--word associations. If sample lag > 0 then the
     * mean value of all sampled statistics for phi[][] is taken.
     * 
     * @return phi multinomial mixture of topic words (K x V)
     */
    public double[][] getPhi() {
        double[][] phi = new double[K][V];
        if (SAMPLE_LAG > 0) {
            for (int k = 0; k < K; k++) {
                for (int w = 0; w < V; w++) {
                    phi[k][w] = phisum[k][w] / numstats;
                }
            }
        } else {
            for (int k = 0; k < K; k++) {
                for (int w = 0; w < V; w++) {
                    phi[k][w] =  (cwt[w][k] + nwk[w][k] + beta) / (cwtsum[k] + nwsum[k] + V * beta);
                }
            }
        }
        return phi;
    }
    
    /**
     * Retrieve estimated topic--word associations. If sample lag > 0 then the
     * mean value of all sampled statistics for phi[][] is taken.
     * 
     * @return phi multinomial mixture of topic words (K x V)
     */
    public int[][] getTestDocTopicAssignmentLattice() {
    	int[][] lat = new int[documents.length*K][];
    	int word, topic;
    	for(int d=0; d<documents.length; d++){
    		for(int k=0; k<K; k++){
	    		lat[d*K + k] = new int[documents[d].length];
	    	}
    		for(int n=0; n<documents[d].length; n++){
    			word = documents[d][n];
    			for(int n2=0; n2<documents[d].length; n2++){
    				if(word == documents[d][n2]){
					topic = z[d][n2];
					lat[d*K + topic][n2]++;
				}
			}
    		}
    	}
        return lat;
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
            System.err.println("\nUSAGE: java LdaGibbsSamplingForInference <testDocumentWordIndexFile> <vocabFile> <K> <alpha> <beta> <trainedCwtFile> <output-tag>\n");
            System.exit(0);
        }
        
        try
    	{
			CorpusReader cr = new CorpusReader(args[0], args[1]);
			int V = cr.getV();
			int M = cr.getD();
			int[][] documents = cr.getDocumentTermIndex();
		
			int K = Integer.parseInt(args[2]);
			double alpha = Float.parseFloat(args[3]);
			double beta = Float.parseFloat(args[4]);
		
			int[][] cwt = ArrayMatrixUtils.readIntMatrixFromFile(args[5]);

			if(V != cwt.length){
				System.out.println("\n Fatal Error: Mismatch between vocab size " + V + " and cwt " + cwt.length);
			}

			System.out.println("Inference using Gibbs Sampling.");
			System.out.println("\nRunning Inference for " + K + " topics for " + M + " documents with " + V + " words vocabulary ...\n");

			LdaGibbsSamplingForInference lda = new LdaGibbsSamplingForInference(documents, V);
			lda.configure(2500, 1000, 100, 10);
			lda.gibbs(K, alpha, beta, cwt);

			double[][] theta = lda.getTheta();

			TopicUtils.saveCSV(theta, args[6]+"-test_theta");
			//TopicUtils.saveCSV(lda.getTestDocTopicAssignmentLattice(), args[6]+"-test_doc_topic-assign-lat");

		}catch (IOException ex){
            System.err.println("\nIOException...\nUSAGE: java LdaGibbsSamplingForInference <testDocumentWordIndexFile> <vocabFile> <K> <alpha> <beta> <trainedCwtFile> <output-tag>\n");
		    System.err.println(ex.getMessage());
		}
    }
}
