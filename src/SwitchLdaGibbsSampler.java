/*
 * SwitchLdaGibbsSampler is free software; you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the Free
 * Software Foundation; either version 2 of the License, or (at your option) any
 * later version.
 *
 * SwitchLdaGibbsSampler is distributed in the hope that it will be useful, but
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
 * Gibbs sampler for estimating the best assignments of topics for words, etities and
 * documents in a corpus. The model and algorithm is introduced in David Newmans' paper
 * "Statistical Entity-Topic Models"
 * (2006).
 * 
 * @author mranahmd
 *
 */
public class SwitchLdaGibbsSampler {

    /**
     * document word data (term lists)
     */
    int[][] documentWords;
    
    /**
     * document entity data (term lists)
     */
    int[][] documentEntities;

    /**
     * vocabulary size of words
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
     * Beta parameter (topic--WordTopicSwitch associations)
     */
    double gamma;

    /**
     * topic assignments for each non-entity.
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
     * cet[i][j] number of instances of entity i (term?) assigned to topic j.
     */
    int[][] cet;

    /**
     * sumcet[j] total number of entities assigned to topic j.
     */
    int[] sumcet;

	/**
     * ctd[i][j] number of words+entities in document i assigned to topic j.
     */
    int[][] ctd;

    /**
     * sumctd[i] total number of words+entities in document i.
     */
    int[] sumctd;
    
	/**
     * nwt[] total number of terms assigned as words by switch/flag in topic i.
     */
    int[] nwt;

	/**
     * net[] total number of terms assigned as entities by switch/flag in topic i.
     */
    int[] net;

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
     * cumulative statistics of psi
     */
    double[] psisum;

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
     *            total vocabulary size
     * @param W
     *            entity vocabulary size
     * @param data
     */
    public SwitchLdaGibbsSampler(int[][] documentWords, int [][] documentEntities, int W, int E) {
        this.documentWords = documentWords;
        this.documentEntities = documentEntities;
        this.E = E;
        this.W = W;
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
        nwt = new int[K];
        net = new int[K];

        // The z_i are are initialised to values in [1,K] to determine the
        // initial state of the Markov chain.

		// Note: This following is a simpler code implementation than maintaining the flag 'x' as mentioned in the paper.
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
                nwt[topic]++;
            }
            N = documentEntities[m].length;
            ze[m] = new int[N];
            for (int n = 0; n < N; n++) {
                int topic = (int) (Math.random() * K);
                ze[m][n] = topic;
                // number of instances of entity i assigned to topic j
                cet[documentEntities[m][n]][topic]++;
                // number of entities in document i assigned to topic j.
                ctd[m][topic]++;
                // total number of entities assigned to topic j.
                sumcet[topic]++;
                net[topic]++;
            }
            // total number of words in document i
            sumctd[m] = documentWords[m].length + documentEntities[m].length;
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
     * @param gamma
     *            symmetric prior parameter on topic--switch associations
     */
    private void gibbs(int K, double alpha, double betaw, double betae, double gamma) {
        this.K = K;
        this.alpha = alpha;
        this.betaw = betaw;
        this.betae = betae;
		this.gamma = gamma;

        // init sampler statistics
        if (SAMPLE_LAG > 0) {
            thetasum = new double[documentWords.length][K];
            phiwsum = new double[K][W];
            phiesum = new double[K][E];
            psisum = new double[K];
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
     *	x_i = 0, therefore
     * 	Sample a topic z_i from the full conditional distribution: p(z_i = j |
     * 	z_-i, w, x_i=0) = (cwt_-i,j(w_i) + betaw)/(sumcwt_-i,j(.) + W * betaw) * 
     *	(ctd_-i,j(d_i) + alpha)/(sumctd_-i,.(d_i) + K * alpha) *
     *	(nwt_-i + gamma)/(nwt_-i + net + 2*gamma)
     * 
     * @param m
     *            document
     * @param n
     *            word
     */
    private int sampleWFullConditional(int m, int n) {

        // remove z_i from the count variables
        int topic = zw[m][n];
        cwt[documentWords[m][n]][topic]--;
        ctd[m][topic]--;
        sumcwt[topic]--;
        sumctd[m]--;
        nwt[topic]--;

        // do multinomial sampling via cumulative method:
        double[] p = new double[K];
        for (int k = 0; k < K; k++) {
            p[k] = (cwt[documentWords[m][n]][k] + betaw) / (sumcwt[k] + W * betaw)
                 * (ctd[m][k] + alpha) / (sumctd[m] + K * alpha)
                 * (nwt[k] + gamma) / (nwt[k] + net[k] + 2*gamma);
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
        cwt[documentWords[m][n]][topic]++;
        ctd[m][topic]++;
        sumcwt[topic]++;
        sumctd[m]++;
        nwt[topic]++;

        return topic;
    }
    
    /**
     *  x_i = 1, therefore 
     *	Sample a topic z_i from the full conditional distribution: p(z_i = j |
     * 	z_-i, w, x_i=1) = (cet_-i,j(w_i) + betae)/(sumcet_-i,j(.) + E * betae) * 
     *	(ctd_-i,j(d_i) + alpha)/(sumcet_-i,.(d_i) + K * alpha) *
     *	(net_-i + gamma)/(nwt + net_-i + 2*gamma)
     * 
     * @param m
     *            document
     * @param n
     *            entity
     */
    private int sampleEFullConditional(int m, int n) {

        // remove z_i from the count variables
        int topic = ze[m][n];
        cet[documentEntities[m][n]][topic]--;
        ctd[m][topic]--;
        sumcet[topic]--;
        sumctd[m]--;
        net[topic]--;

        // do multinomial sampling via cumulative method:
        double[] p = new double[K];
        for (int k = 0; k < K; k++) {
            p[k] = (cet[documentEntities[m][n]][k] + betae) / (sumcet[k] + E * betae)
                 * (ctd[m][k] + alpha) / (sumctd[m] + K * alpha)
                 * (net[k] + gamma) / (net[k] + nwt[k] + 2*gamma);
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
        cet[documentEntities[m][n]][topic]++;
        ctd[m][topic]++;
        sumcet[topic]++;
        sumctd[m]++;
        net[topic]++;

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
        for (int k = 0; k < K; k++) {
            psisum[k] += (net[k] + gamma) / (net[k] + nwt[k] + 2*gamma);
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
     * Retrieve estimated topic--switch associations. If sample lag > 0 then the
     * mean value of all sampled statistics for psi[] is taken.
     * 
     * @return phie multinomial mixture of topic words (K x E)
     */
    public double[] getPsi() {
        double[] psi = new double[K];
        if (SAMPLE_LAG > 0) {
            for (int k = 0; k < K; k++) {
                psi[k] = psisum[k] / numstats;
            }
        } else {
            for (int k = 0; k < K; k++) {
                psi[k] = (net[k] + gamma) / (net[k] + nwt[k] + 2*gamma);
            }
        }
        return psi;
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
     * 
     * @return net (K)
     */
    public int[] getNwt() {
        return nwt;
    }
    
     /**
     * 
     * @return net (K)
     */
    public int[] getNet() {
        return net;
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
		if (args.length != 9)
        {
            System.err.println("\nUSAGE: java SwitchLdaGibbsSampler <documentWordIndexFile> <wordVocabFile> <documentEntityIndexFile> <entityVocabFile> <K> <alpha> <betaw> <betae> <gamma>\n");
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
			double gamma = Float.parseFloat(args[8]);

			System.out.println("Switch Latent Dirichlet Allocation using Gibbs Sampling.");
			System.out.println("\nTraining " + K + " SwitchLDA topics for " + M + " documents with " + W + " words vocabulary and " + E + " entities vocabulary ...\n");

			SwitchLdaGibbsSampler swlda = new SwitchLdaGibbsSampler(documentWords, documentEntities, W, E);
			swlda.configure(2500, 1000, 100, 10);
			swlda.gibbs(K, alpha, betaw, betae, gamma);

			double[][] theta = swlda.getTheta();
			double[][] phiw = swlda.getPhiW();
			double[][] phie = swlda.getPhiE();
			double[] psi = swlda.getPsi();
			int[][] cwt = swlda.getCwt();
			int[][] cet = swlda.getCet();
			int[] nwt = swlda.getNwt();
			int[] net = swlda.getNet();
			
			TopicUtils.saveCSV(theta,"./SwitchLDA_theta");
			TopicUtils.saveCSV(phiw,"./SwitchLDA_phiw");
			TopicUtils.saveCSV(phie,"./SwitchLDA_phie");
			TopicUtils.saveCSV(psi,"./SwitchLDA_psi");
			TopicUtils.saveCSV(cwt,"./SwitchLDA_cwt");
			TopicUtils.saveCSV(cet,"./SwitchLDA_cet");
			TopicUtils.saveCSV(nwt,"./SwitchLDA_nwt");
			TopicUtils.saveCSV(net,"./SwitchLDA_net");
			TopicUtils.saveTopTopicWords(phiw, wcr.getVocab(), "./SwitchLDA_top-topic-words", 10);
			TopicUtils.saveTopTopicWords(phie, ecr.getVocab(), "./SwitchLDA_top-topic-entities", 10);
		}catch (IOException ex){
			System.err.println("\nIOException....\nUSAGE: java SwitchLdaGibbsSampler <documentWordIndexFile> <wordVocabFile> <documentEntityIndexFile> <entityVocabFile> <K> <alpha> <betaw> <betae> <gamma>\n");
			System.err.println(ex.getMessage());
		}
	}
}
