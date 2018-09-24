# viralrank_centrality
ViralRank is a centrality measure that identifies the most influential nodes in complex networks, defined in 
[F. Iannelli, M.S. Mariani, I.M. Sokolov, "Influencers identification in complex networks through reaction-diffusion dynamics"] (https://arxiv.org/abs/1803.01212).

Viralrank is defined as (note that the minus sign implies negative scores): 

v(lambda) = - \sum_j (D_{ij} + D_{ji})/N,
        
where D_{ij} is the random-walk effective distance [F. Iannelli at al. Phys. Rev. E 95, 012313 (2017)] from i to j and lambda is a parameter that plays the role of the inverse temperature (see next). The latter is defined as:

lambda = log[ (beta-mu)/alpha ] - gamma

where beta, mu and alpha are the infection, recovery and diffusion rates, respectively, while gamma is the Euler constant.
For contact networks its value is set equal to (approximately) zero corresponding to a high temperature expansion, while it can be specifically tuned for metapopulation networks uisng the above definition.
                  
# required modules and basic usage

---------
Required modules:

Python:   tested for: 2.7.13.  
numpy:    tested for: 1.11.0.  
networkx: tested for: 1.14.5.   

After importing the module with:

import vrank as vr

ViralRank scores can be computed for the networkx graph 'g' as (for contact networks set e.g. lambda = 0.0001):

v = vr.ViralRank(g).value(lambda) 

The corresponding nodes' ranking is obtained as:

v = vr.ViralRank(g).rank(lambda) 
