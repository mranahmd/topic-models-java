/*
 * CorrLda1GibbsSampler is free software; you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the Free
 * Software Foundation; either version 2 of the License, or (at your option) any
 * later version.
 *
 * CorrLda1GibbsSampler is distributed in the hope that it will be useful, but
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
 * Gibbs sampler for estimating the best assignments of topics for words, entities and
 * documents in a corpus. The model and algorithm is introduced in David Newmans' paper
 * "Statistical Entity-Topic Models"
 * (2006).
 * 
 * @author mranahmd
 *
 */
public class CorrLda1GibbsSampler {

    /**
     * document data (term lists)
     */
    int[][] documentWords;
    
    /**
     * document data (entity lists)
     */
    int[][] documentEntities;

    /**
     * vocabulary size of words other than entities
     */
    int W;
    
    /**
     * vocabulary size of entities
     */
    int E;

    /**
     * number of topics
     */
    int K;

    /**
     * Dirichlet parameter (document--topic associations)
     */
    double alpha;

    /**
     * Dirichlet parameter (topic--word associations)
     */
    double betaw;
    
    /**
     * Dirichlet parameter (topic--entity associations)
     */
    double betae;

    /**
     * topic assignments for each word.
     */
    int zw[][];

    /**
     * topic assignments for each entity.
     */
    int ze[][];

    /**
     * cwt[i][j] number of instances of word i (term?) assigned to topic j.
     */
    int[][] cwt;

    /**
     * sumcwt[j] total number of words assigned to topic j.
     */
    int[] sumcwt;
    
    /**
     * cet[i][j] number of instances of entity i assigned to topic j.
     */
    int[][] cet;

    /**
     * sumcet[j] total number of entities assigned to topic j.
     */
    int[] sumcet;

	/**
     * ctd[i][j] number of terms (words+entities) in document i assigned to topic j.
     */
    int[][] ctd;

    /**
     * sumctd[i] total number of words+entities in document i.
     */
    int[] sumctd;

    /**
     * cumulative statistics of theta
     */
    double[][] thetasum;

    /**
     * cumulative statistics of phiw
     */
    double[][] phiwsum;
    
    /**
     * cumulative statistics of phie
     */
    double[][] phiesum;

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
    public CorrLda1GibbsSampler(int[][] documentWords, int[][] documentEntities, int W, int E) {
        this.documentWords = documentWords;
        this.documentEntities = documentEntities;
        this.W = W;
        this.E = E;		
    }

    /**
     * Initialisation: Must start with an assignment of observations to topics ?
     * Many alternatives are possible, I chose to perform random assignments
     * with equal probabilities
     * 
     * @param K
     *            number of topics
     * @return zw assignment of topics to words
     *		   ze assignment of topics to entities
     */
    public void initialState(int K) {
        int i;

        int M = documentWords.length;

        // initialise count variables.
        cwt = new int[W][K];
        cet = new int[E][K];
        ctd = new int[M][K];
        sumcwt = new int[K];
        sumcet = new int[K];
        sumctd = new int[M];

        // The z_i are are initialised to values in [1,K] to determine the
        // initial state of the Markov chain.

        zw = new int[M][];
        ze = new int[M][];
        for (int m = 0; m < M; m++) {
            int N = documentWords[m].length;
            zw[m] = new int[N];
            for (int n = 0; n < N; n++) {
                int topic = (int) (Math.random() * K);
                zw[m][n] = topic;
                // number of instances of word i assigned to topic j
                cwt[documentWords[m][n]][topic]++;
                // number of words in document i assigned to topic j.
                ctd[m][topic]++;
                // total number of words assigned to topic j.
                sumcwt[topic]++;
            }
            N = documentEntities[m].length;
            ze[m] = new int[N];
            for (int n = 0; n < N; n++) {
                int topic = (int) (Math.random() * K);
                ze[m][n] = topic;
                // number of instances of entity i assigned to topic j
                cet[documentEntities[m][n]][topic]++;
                // total number of entities assigned to topic j.
                sumcet[topic]++;
            }
            // total number of words in document i
            sumctd[m] = documentWords[m].length;
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
     * @param betaw
     *            symmetric prior parameter on topic--word associations
     * @param betae
     *            symmetric prior parameter on topic--entity associations
     */
    private void gibbs(int K, double alpha, double betaw, double betae) {
        this.K = K;
        this.alpha = alpha;
        this.betaw = betaw;
        this.betae = betae;

        // init sampler statistics
        if (SAMPLE_LAG > 0) {
            thetasum = new double[documentWords.length][K];
            phiwsum = new double[K][W];
            phiesum = new double[K][E];
            numstats = 0;
        }

        // initial state of the Markov chain:
        initialState(K);

        System.out.println("Sampling " + ITERATIONS
            + " iterations with burn-in of " + BURN_IN + " (B/S="
            + THIN_INTERVAL + ").");

        for (int i = 0; i < ITERATIONS; i++) {

            // for all z_i
            for (int m = 0; m < zw.length; m++) {
                for (int n = 0; n < zw[m].length; n++) {

                    // (z_i = z[m][n])
                    // sample from p(z_i|z_-i, w)
                    int topic = sampleWFullConditional(m, n);
                    zw[m][n] = topic;
                }
                for (int n = 0; n < ze[m].length; n++) {

                    // (z_i = z[m][n])
                    // sample from p(z_i|z_-i, w)
                    int topic = sampleEFullConditional(m, n);
                    ze[m][n] = topic;
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
     * Sample a topic zw_i from the full conditional distribution: p(zw_i = j |
     * zw_-i, w) = (cwt_-i,j(w_i) + betaw)/(sumcwt_-i,j(.) + W * betaw) * (ctd_-i,j(d_i) +
     * alpha)/(sumctd_-i,.(d_i) + K * alpha)
     * 
     * @param m
     *            document
     * @param n
     *            word
     */
    private int sampleWFullConditional(int m, int n) {

        // remove zw_i from the count variables
        int topic = zw[m][n];
        cwt[documentWords[m][n]][topic]--;
        ctd[m][topic]--;
        sumcwt[topic]--;
        sumctd[m]--;

        // do multinomial sampling via cumulative method:
        double[] p = new double[K];
        for (int k = 0; k < K; k++) {
            p[k] = (cwt[documentWords[m][n]][k] + betaw) / (sumcwt[k] + W * betaw)
                * (ctd[m][k] + alpha) / (sumctd[m] + K * alpha);
        }
        // cumulate multinomial parameters
        for (int k = 1; k < p.length; k++) {
            p[k] += p[k - 1];
        }
        // scaled sample because of unnormalised p[]
        double u = Math.random() * p[p.length - 1];
        for (topic = 0; topic < p.length; topic++) {
            if (u < p[topic])
                break;
        }

        // add newly estimated zw_i to count variables
        cwt[documentWords[m][n]][topic]++;
        ctd[m][topic]++;
        sumcwt[topic]++;
        sumctd[m]++;

        return topic;
    }
    
    /**
     * Sample a topic ze_i from the full conditional distribution: p(ze_i = j |
     * ze, e, zw_-i) = (cet_-i,j(w_i) + betae)/(sumcet_-i,j(.) + E * betae) * (ctd_i,j(d_i))/(W)
     * 
     * @param m
     *            document
     * @param n
     *            entity
     */
    private int sampleEFullConditional(int m, int n) {

        // remove ze_i from the count variables
        int topic = ze[m][n];
        cet[documentEntities[m][n]][topic]--;
        sumcet[topic]--;

        // do multinomial sampling via cumulative method:
        double[] p = new double[K];
        for (int k = 0; k < K; k++) {
            p[k] = (cet[documentEntities[m][n]][k] + betae) / (sumcet[k] + E * betae)
                * (ctd[m][k]) / (documentWords[m].length);
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

        // add newly estimated ze_i to count variables
        cet[documentEntities[m][n]][topic]++;
        sumcet[topic]++;
        
        return topic;
    }

    /**
     * Add to the statistics the values of theta and phi for the current state.
     */
    private void updateParams() {
        for (int m = 0; m < documentWords.length; m++) {
            for (int k = 0; k < K; k++) {
                thetasum[m][k] += (ctd[m][k] + alpha) / (sumctd[m] + K * alpha);
            }
        }
        for (int k = 0; k < K; k++) {
            for (int w = 0; w < W; w++) {
                phiwsum[k][w] += (cwt[w][k] + betaw) / (sumcwt[k] + W * betaw);
            }
        }
        for (int k = 0; k < K; k++) {
            for (int e = 0; e < E; e++) {
                phiesum[k][e] += (cet[e][k] + betae) / (sumcet[k] + E * betae);
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
        double[][] theta = new double[documentWords.length][K];

        if (SAMPLE_LAG > 0) {
            for (int m = 0; m < documentWords.length; m++) {
                for (int k = 0; k < K; k++) {
                    theta[m][k] = thetasum[m][k] / numstats;
                }
            }

        } else {
            for (int m = 0; m < documentWords.length; m++) {
                for (int k = 0; k < K; k++) {
                    theta[m][k] = (ctd[m][k] + alpha) / (sumctd[m] + K * alpha);
                }
            }
        }

        return theta;
    }

    /**
     * Retrieve estimated topic--word associations. If sample lag > 0 then the
     * mean value of all sampled statistics for phiw[][] is taken.
     * 
     * @return phiw multinomial mixture of topic words (K x W)
     */
    public double[][] getPhiW() {
        double[][] phiw = new double[K][W];
        if (SAMPLE_LAG > 0) {
            for (int k = 0; k < K; k++) {
                for (int w = 0; w < W; w++) {
                    phiw[k][w] = phiwsum[k][w] / numstats;
                }
            }
        } else {
            for (int k = 0; k < K; k++) {
                for (int w = 0; w < W; w++) {
                    phiw[k][w] = (cwt[w][k] + betaw) / (sumcwt[k] + W * betaw);
                }
            }
        }
        return phiw;
    }
    
    /**
     * Retrieve estimated topic--entity associations. If sample lag > 0 then the
     * mean value of all sampled statistics for phie[][] is taken.
     * 
     * @return phie multinomial mixture of topic words (K x E)
     */
    public double[][] getPhiE() {
        double[][] phie = new double[K][E];
        if (SAMPLE_LAG > 0) {
            for (int k = 0; k < K; k++) {
                for (int e = 0; e < E; e++) {
                    phie[k][e] = phiesum[k][e] / numstats;
                }
            }
        } else {
            for (int k = 0; k < K; k++) {
                for (int e = 0; e < E; e++) {
                    phie[k][e] = (cet[e][k] + betae) / (sumcet[k] + E * betae);
                }
            }
        }
        return phie;
    }
    
     /**
     * 
     * @return cwt (K x V)
     */
    public int[][] getCwt() {
        return cwt;
    }
    
     /**
     * 
     * @return cet (K x V)
     */
    public int[][] getCet() {
        return cet;
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
		if (args.length != 8)
        {
            System.err.println("\nUSAGE: java CorrLda1GibbsSampler <documentWordIndexFile> <wordVocabFile> <documentEntityIndexFile> <entityVocabFile> <K> <alpha> <betaw> <betae>\n");
            System.exit(0);
        }
        
        try
    	{
			CorpusReader wcr = new CorpusReader(args[0], args[1]);
			int W = wcr.getV();
			int M = wcr.getD();
			int[][] documentWords = wcr.getDocumentTermIndex();
		
			CorpusReader ecr = new CorpusReader(args[2], args[3]);
			int E = ecr.getV();
			int[][] documentEntities = ecr.getDocumentTermIndex();
			if(M != ecr.getD()){
				System.out.println("\n Fatal Error: Mismatch between word documents size " + M + " and entity document size " + ecr.getD());
			}
		
			int K = Integer.parseInt(args[4]);
			double alpha = Float.parseFloat(args[5]);
			double betaw = Float.parseFloat(args[6]);
			double betae = Float.parseFloat(args[7]);

			System.out.println("Corr Latent Dirichlet Allocation 1 using Gibbs Sampling.");
			System.out.println("\nTraining " + K + " CorrLDA1 topics for " + M + " documents with " + W + " words vocabulary and " + E + " entities vocabulary ...\n");

			CorrLda1GibbsSampler corrlda1 = new CorrLda1GibbsSampler(documentWords, documentEntities, W, E);
			corrlda1.configure(2500, 1000, 100, 10);
			corrlda1.gibbs(K, alpha, betaw, betae);

			double[][] theta = corrlda1.getTheta();
			double[][] phiw = corrlda1.getPhiW();
			double[][] phie = corrlda1.getPhiE();
			int[][] cwt = corrlda1.getCwt();
			int[][] cet = corrlda1.getCet();

			TopicUtils.saveCSV(theta,"./CorrLDA1_theta");
			TopicUtils.saveCSV(phiw,"./CorrLDA1_phiw");
			TopicUtils.saveCSV(phie,"./CorrLDA1_phie");
			TopicUtils.saveCSV(cwt,"./CorrLDA1_cwt");
			TopicUtils.saveCSV(cet,"./CorrLDA1_cet");
			TopicUtils.saveTopTopicWords(phiw, wcr.getVocab(), "./CorrLDA1_top-topic-words", 10);
			TopicUtils.saveTopTopicWords(phie, ecr.getVocab(), "./CorrLDA1_top-topic-entities", 10);

		}catch (IOException ex){
			System.err.println("\nIOException....\nUSAGE: java CorrLda1GibbsSampler <documentWordIndexFile> <wordVocabFile> <documentEntityIndexFile> <entityVocabFile> <K> <alpha> <betaw> <betae>\n");
			System.err.println(ex.getMessage());
		}
    }
}
