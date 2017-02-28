# NetICS
Network-based integration of multi-omics data for prioritizing cancer genes

NetICS performs a per sample bidirectional network diffusion technique to prioritize genes based on their proximity to genetically aberrant and differentially expressed genes. It provides rank aggregation techniques for integrating the sample-specific gene lists into an overall ranked list of genes.

The method should be called as follows:

```
ranked_list_genes = netics_fun( 'mutation_data_breast.txt', F, F_opp, 'RANK_AGGREG=SUM', 'UP_DOWN=3', 'breastDEgenes.txt', 'network_genes.txt');
```

'mutation_data_breast.txt' --> a file with two columns that maps every gene (1^st column) to the samples that it is genetically aberrant (2<sup>nd</sup> column).

'breastDEgenes.txt' --> a file containing a list of gene names detected as differentially expressed.

'network_genes.txt' --> a file containing the list of network gene names. They should be in the same order as in the diffused matrices F and F_opp.

Example files for 'breastDEgenes.txt' and 'network_genes.txt' are given.

RANK_AGGREG can take values "SUM", "MEDIAN" or "RRA". "SUM" computes the summation of the per sample ranks and "MEDIAN" computes the median. "RRA" computes the robust rank aggregation technique as described in (Kolde et al, 2012). The matlab code for implementing the RRA method can be found at http://ch.mathworks.com/matlabcentral/fileexchange/41835-rank-aggregation. You will need to include the files betaScores.m, correctBetaPvalues.m, rhoScores.m and thresholdBetaScore.m.

UP_DOWN can take values from 1 to 3. "1" performs diffusion from the genetically aberrant genes, "2" performs diffusion from the differentially expressed genes, and "3" performs bidirectional network diffusion (from both aberrant and differentially expressed genes).

F and F_opp are precomputed diffused matrices on a given network.



### Contributions
- [Christos Dimitrakopoulos](https://www.bsse.ethz.ch/cbg/group/people/person-detail.html?persid=197642)
- [Niko Beerenwinkel](http://www.bsse.ethz.ch/cbg/group/people/person-detail.html?persid=149417)


###Contact
```
Christos Dimitrakopoulos
christos.dimitrakopoulos (at) bsse.ethz.ch
```
