/*
 * Lda.java is free software; you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the Free
 * Software Foundation; either version 2 of the License, or (at your option) any
 * later version.
 *
 * Lda.java is distributed in the hope that it will be useful, but
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
 * 	Multispeech - LORIA
 * 	Nancy, France
 *
 *  Acknowledgement: Gregor Heinrich (gregor :: arbylon : net) 
 *  (This file is uses the Java code of Gregor Heinrich (gregor@arbylon.net)
 *    http://www.arbylon.net/projects/LdaGibbsSampler.java
 *  part of the org.knowceans experimental software packages.)
 * 
 */

import java.io.*;
import org.knowceans.gibbstest.LdaGibbsSampler;
import java.lang.reflect.Method;
import java.lang.reflect.InvocationTargetException;
import java.lang.reflect.Field;

/**
 * Class to test LDA topic model in org.knowceans.gibbstest package.
 * 
 * @author mranahmd
 *
 */
public class TestLda { 
	public static void main(String[] args) {
	if (args.length != 5)
        {
            System.err.println("\nUSAGE: java LdaGibbsSampler <documentTermIndexFile> <vocabFile> <K> <alpha> <beta>\n");
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

			System.out.println("Latent Dirichlet Allocation using Gibbs Sampling.");
			System.out.println("\nTraining " + K + " LDA topics for " + M + " documents with " + V + " words vocabulary ...\n");

			LdaGibbsSampler lda = new LdaGibbsSampler(documents, V);
			Method ldaGibbs = LdaGibbsSampler.class.getDeclaredMethod("gibbs");
			ldaGibbs.setAccessible(true);
			Field ldaNw = LdaGibbsSampler.class.getDeclaredField("nw"); 
			ldaNw.setAccessible(true);
			
			lda.configure(2500, 1000, 100, 10);
			ldaGibbs.invoke(lda, K, alpha, beta);

			double[][] theta = lda.getTheta();
			double[][] phi = lda.getPhi();
			int[][] cwt = (int[][]) ldaNw.get(lda);

			TopicUtils.saveCSV(theta,"./LDA_theta");
			TopicUtils.saveCSV(phi,"./LDA_phi");
			TopicUtils.saveCSV(cwt,"./LDA_cwt");
			TopicUtils.saveTopTopicWords(phi, cr.getVocab(), "./LDA_top-topic-words", 10);

		}catch (IOException ex){
           		System.err.println("\nIOException ...\nUSAGE: java LdaGibbsSampler <documentTermIndexFile> <vocabFile> <K> <alpha> <beta>\n");
		        System.err.println(ex.getMessage());
		}catch (NoSuchMethodException ex){
		        System.err.println(ex.getMessage());
		}catch (IllegalAccessException ex){
		        System.err.println(ex.getMessage());
		}catch (InvocationTargetException ex){
		        System.err.println(ex.getMessage());
		}catch (NoSuchFieldException ex){
		        System.err.println(ex.getMessage());
		}
    }
}
