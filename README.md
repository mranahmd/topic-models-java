# topic-models-java
This repository has Java implementation of some interesting topic models based on [Latent Dirichlet Allocation (LDA)](https://en.wikipedia.org/wiki/Topic_model) [1].

Currently it includes the implementation of Entity-Topic models [2][3], the Author-Topic model [4] and Special Words with Background (SWB) model [5]. For each of the topic models, the model parameters are estimated using [Gibbs Sampling](https://en.wikipedia.org/wiki/Gibbs_sampling).

The code shared in this repository is inspired from the simple single file implementation of Gibbs Sampler for LDA shared [here](http://www.arbylon.net/resources.html).

For more topic models and their implementations check [here](http://mimno.infosci.cornell.edu/topics.html).

References: <br>
[1] David M. Blei, Andrew Y. Ng, and Michael I. Jordan. 2003. Latent dirichlet allocation. J. Mach. Learn. Res. 3 (March 2003), 993-1022. <br>
[2] David Newman, Chaitanya Chemudugunta, Padhraic Smyth. Statistical entity-topic models. KDD (2006). <br>
[3] Linmei Hu, Juanzi Li, Zhihui Li, Chao Shao, Zhixing Li. Incorporating Entities in News Topic Modeling. NLPCC 2013, CCIS 400. <br>
[4] Michal Rosen-Zvi, Thomas Griffiths, Mark Steyvers, and Padhraic Smyth. 2004. The author-topic model for authors and documents. UAI '04. <br>
[5] Chaitanya Chemudugunta, Padhraic Smyth, Mark Steyvers. Modeling General and Specific Aspects of Documents with a Probabilistic Topic Model. NIPS 2006: 241-248. <br>

