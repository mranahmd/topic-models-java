/*
 * CorrLda2GibbsSampler is free software; you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the Free
 * Software Foundation; either version 2 of the License, or (at your option) any
 * later version.
 *
 * CorrLda2GibbsSampler is distributed in the hope that it will be useful, but
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
public class CorrLda2GibbsSampler {

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
     * number of word topics
     */
    int KW;
    
    /**
     * number of entity topics
     */
    int KE;
    
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
     * Dirichlet parameter (topic--supertopic associations)
     */
    double gamma;

    /**
     * topic assignments for each word.
     */
    int zw[][];

    /**
     * topic assignments for each entity.
     */
    int ze[][];
    
    /**
     * supertopic assignments for each entity.
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
     * cet[i][j] number of instances of entity i (term?) assigned to topic j.
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
     * ctt[i][j] number of topics per supertopic.
     */
    int[][] ctt;
    
    /**
     * sumctt[i] total number of topics per supertopic.
     */
    int[] sumctt;
        
    /**    
    
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
    double[][] psisum;

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
    public CorrLda2GibbsSampler(int[][] documentWords, int[][] documentEntities, int W, int E) {
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
    public void initialState(int KW, int KE) {
        int i;

        int M = documentWords.length;

        // initialise count variables.
        cwt = new int[W][KW];
        cet = new int[E][KE];
        ctd = new int[M][KW];
        ctt = new int[KW][KE];        
        sumcwt = new int[KW];
        sumcet = new int[KE];
        sumctd = new int[M];
        sumctt = new int[KW];

        // The z_i are are initialised to values in [1,K] to determine the
        // initial state of the Markov chain.

        zw = new int[M][];
        ze = new int[M][];
        x = new int[M][];
        for (int m = 0; m < M; m++) {
            int N = documentWords[m].length;
            zw[m] = new int[N];
            for (int n = 0; n < N; n++) {
                int topic = (int) (Math.random() * KW);
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
            x[m] = new int[N];
            for (int n = 0; n < N; n++) {                
                int topicze = (int) (Math.random() * KE);
                ze[m][n] = topicze;
                // number of instances of entity i assigned to topic j
                cet[documentEntities[m][n]][topicze]++;

                // total number of entities assigned to topic j.
                sumcet[topicze]++;
                
                int supertopic = (int) (Math.random() * KW);
                x[m][n] = supertopic;
                ctt[supertopic][topicze]++;
                sumctt[supertopic]++;
            
            }
            // total number of words in document i
            sumctd[m] = documentWords[m].length;// + documentEntities[m].length;
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
    private void gibbs(int KW, int KE, double alpha, double betaw, double betae, double gamma) {
        this.KW = KW;
        this.KE = KE;
        this.alpha = alpha;
        this.betaw = betaw;
        this.betae = betae;
        this.gamma = gamma;

        // init sampler statistics
        if (SAMPLE_LAG > 0) {
            thetasum = new double[documentWords.length][KW];
            phiwsum = new double[KW][W];
            phiesum = new double[KE][E];
            psisum = new double[KW][KE];
            numstats = 0;
        }

        // initial state of the Markov chain:
        initialState(KW, KE);

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
                    int topics[] = sampleEFullConditional(m, n);
                    x[m][n] = topics[0];	//supertopic
                    ze[m][n] = topics[1];	//topic
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
     * alpha)/(sumctd_-i,.(d_i) + KW * alpha)
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
        double[] p = new double[KW];
        for (int k = 0; k < KW; k++) {
            p[k] = (cwt[documentWords[m][n]][k] + betaw) / (sumcwt[k] + W * betaw)
                * (ctd[m][k] + alpha) / (sumctd[m] + KW * alpha);
        }
        // cumulate multinomial parameters
        for (int k = 1; k < p.length; k++) {
            p[k] += p[k - 1];
        }
        // scaled sample because of unnormalised p[]
        double u = Math.random() * p[KW - 1];
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
     * ze, e, zw_-i) = (cet_-i,j(e_i) + betae)/(sumcet_-i,j(.) + E * betae) * (ctd_-i,j(d_i))/(W) * (ctt_-st,t + gamma)/(sumctt- + KE * gamma)
     * 
     * @param m
     *            document
     * @param n
     *            entity
     */
    private int[] sampleEFullConditional(int m, int n) {

        // remove ze_i from the count variables
        int topic = ze[m][n];
        cet[documentEntities[m][n]][topic]--;
        sumcet[topic]--;
        
        // remove x_i from the count variables
        int supertopic = x[m][n];
        ctt[supertopic][topic]--;
        sumctt[supertopic]--;
        
        // do multinomial sampling via cumulative method:
        double[][] p = new double[KW][KE];
        for (int xk = 0; xk < KW; xk++) {
        	for (int zk = 0; zk < KE; zk++) {
            	p[xk][zk] = (cet[documentEntities[m][n]][zk] + betae) / (sumcet[zk] + E * betae)
                			* (ctd[m][xk]) / (documentWords[m].length)
                			* (ctt[xk][zk] + gamma) / (sumctt[xk] + KE * gamma);
            }
        }
        // cumulate multinomial parameters
        double cump = 0;
        for (int xk = 0; xk < KW; xk++) {
        	for (int zk = 0; zk < KE; zk++) {            
        		p[xk][zk] += cump;
        		cump = p[xk][zk];
        	}
        }
        // scaled sample because of unnormalised p[]
        double u = Math.random() * p[KW-1][KE-1];
        boolean breakFlag = false;
        for (int xk = 0; xk < KW; xk++) {
        	for (int zk = 0; zk < KE; zk++) {              
        		if (u < p[xk][zk]){
        			supertopic = xk;
        			topic = zk;
		            breakFlag = true;
		            break;
		        }
            }
            if (breakFlag){
	            break;
		    }
        }

        // add newly estimated ze_i to count variables
        cet[documentEntities[m][n]][topic]++;
        sumcet[topic]++;        
        ctt[supertopic][topic]++;
        sumctt[supertopic]++;

   		return new int[] {supertopic, topic};
    }

    /**
     * Add to the statistics the values of theta and phi for the current state.
     */
    private void updateParams() {
        for (int m = 0; m < documentWords.length; m++) {
            for (int k = 0; k < KW; k++) {
                thetasum[m][k] += (ctd[m][k] + alpha) / (sumctd[m] + KW * alpha);
            }
        }
        for (int k = 0; k < KW; k++) {
            for (int w = 0; w < W; w++) {
                phiwsum[k][w] += (cwt[w][k] + betaw) / (sumcwt[k] + W * betaw);
            }
        }
        for (int k = 0; k < KE; k++) {
            for (int e = 0; e < E; e++) {
                phiesum[k][e] += (cet[e][k] + betae) / (sumcet[k] + E * betae);
            }
        }
        for (int supertopic = 0; supertopic < KW; supertopic++) {
            for (int topic = 0; topic < KE; topic++) {
                psisum[supertopic][topic] += (ctt[supertopic][topic] + gamma) / (sumctt[supertopic] + KE * gamma);
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
        double[][] theta = new double[documentWords.length][KW];

        if (SAMPLE_LAG > 0) {
            for (int m = 0; m < documentWords.length; m++) {
                for (int k = 0; k < KW; k++) {
                    theta[m][k] = thetasum[m][k] / numstats;
                }
            }

        } else {
            for (int m = 0; m < documentWords.length; m++) {
                for (int k = 0; k < KW; k++) {
                    theta[m][k] = (ctd[m][k] + alpha) / (sumctd[m] + KW * alpha);
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
        double[][] phiw = new double[KW][W];
        if (SAMPLE_LAG > 0) {
            for (int k = 0; k < KW; k++) {
                for (int w = 0; w < W; w++) {
                    phiw[k][w] = phiwsum[k][w] / numstats;
                }
            }
        } else {
            for (int k = 0; k < KW; k++) {
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
        double[][] phie = new double[KE][E];
        if (SAMPLE_LAG > 0) {
            for (int k = 0; k < KE; k++) {
                for (int e = 0; e < E; e++) {
                    phie[k][e] = phiesum[k][e] / numstats;
                }
            }
        } else {
            for (int k = 0; k < KE; k++) {
                for (int e = 0; e < E; e++) {
                    phie[k][e] = (cet[e][k] + betae) / (sumcet[k] + E * betae);
                }
            }
        }
        return phie;
    }
    
    /**
     * Retrieve estimated supertopic--topic associations. If sample lag > 0 then the
     * mean value of all sampled statistics for psi[][] is taken.
     * 
     * @return psi multinomial mixture of topic words (KW x KE)
     */
    public double[][] getPsi() {
        double[][] psi = new double[KW][KE];
        if (SAMPLE_LAG > 0) {
     		for (int supertopic = 0; supertopic < KW; supertopic++) {
            	for (int topic = 0; topic < KE; topic++) {
                    psi[supertopic][topic] = psisum[supertopic][topic] / numstats;
                }
            }
        } else {
     		for (int supertopic = 0; supertopic < KW; supertopic++) {
            	for (int topic = 0; topic < KE; topic++) {
                	psi[supertopic][topic] = (ctt[supertopic][topic] + gamma) / (sumctt[supertopic] + KE * gamma);
            	}
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
     * @return cet (K x V)
     */
    public int[][] getCtt() {
        return ctt;
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
	if (args.length != 10)
        {
            System.err.println("\nUSAGE: java CorrLda2GibbsSampler <documentWordIndexFile> <wordVocabFile> <documentEntityIndexFile> <entityVocabFile> <KW> <KE> <alpha> <betaw> <betae> <gamma>\n");
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

			int KW = Integer.parseInt(args[4]);
			int KE = Integer.parseInt(args[5]);
			double alpha = Float.parseFloat(args[6]);
			double betaw = Float.parseFloat(args[7]);
			double betae = Float.parseFloat(args[8]);
			double gamma = Float.parseFloat(args[9]);

			System.out.println("Corr Latent Dirichlet Allocation 2 using Gibbs Sampling.");
			System.out.println("\nTraining (" + KW + " word | " + KE + " entity) CorrLDA2 topics for " + M + " documents with " + W + " words vocabulary and " + E + " entities vocabulary ...\n");

			CorrLda2GibbsSampler corrlda2 = new CorrLda2GibbsSampler(documentWords, documentEntities, W, E);
			corrlda2.configure(2500, 1000, 100, 10);
			corrlda2.gibbs(KW, KE, alpha, betaw, betae, gamma);

			double[][] theta = corrlda2.getTheta();
			double[][] phiw = corrlda2.getPhiW();
			double[][] phie = corrlda2.getPhiE();
			double[][] psi = corrlda2.getPsi();
			int[][] cwt = corrlda2.getCwt();
			int[][] cet = corrlda2.getCet();
			int[][] ctt = corrlda2.getCtt();

			TopicUtils.saveCSV(theta,"./CorrLDA2_theta");
			TopicUtils.saveCSV(phiw,"./CorrLDA2_phiw");
			TopicUtils.saveCSV(phie,"./CorrLDA2_phie");
			TopicUtils.saveCSV(psi,"./CorrLDA2_psi");
			TopicUtils.saveCSV(cwt,"./CorrLDA2_cwt");
			TopicUtils.saveCSV(cet,"./CorrLDA2_cet");
			TopicUtils.saveCSV(ctt,"./CorrLDA2_ctt");
			TopicUtils.saveTopTopicWords(phiw, wcr.getVocab(), "./CorrLDA2_top-topic-words", 10);
			TopicUtils.saveTopTopicWords(phie, ecr.getVocab(), "./CorrLDA2_top-topic-entities", 10);
			TopicUtils.saveTopTopicsOfSupertopics(psi, "./CorrLDA2_top-Supertopic-topics", 5);

		}catch (IOException ex){
			System.err.println("\nIOException....\nUSAGE: java CorrLda2GibbsSampler <documentWordIndexFile> <wordVocabFile> <documentEntityIndexFile> <entityVocabFile> <K> <alpha> <betaw> <betae>\n");
			System.err.println(ex.getMessage());
		}
    }
}
