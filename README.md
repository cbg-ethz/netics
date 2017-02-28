# NetICS
Network-based integration of multi-omics data for prioritizing cancer genes

NetICS performs a per sample bidirectional network diffusion technique to prioritize genes based on their proximity to genetically aberrant and differentially expressed genes. It provides rank aggregation techniques for integrating the sample-specific gene lists into an overall ranked list of genes.

The method should be called as follows:

```
>> ranked_list_genes = netics_fun( 'mutation_data_breast.txt', F, F_opp, 'RANK_AGGREG=SUM', 'UP_DOWN=3', 'breastDEgenes.txt', 'network_genes.txt');
```

'mutation_data_breast.txt' --> input file that contains the genetically aberrant genes of each sample. It contains two columns that map every gene (1<sup>st</sup> column) to the samples that it is genetically aberrant (2<sup>nd</sup> column).

'breastDEgenes.txt' --> input file that contains a list of gene names detected as differentially expressed between the tumor and normal samples.

'network_genes.txt' --> input file that contains the list of network gene names. They should be in the same order as in the rows of the diffused matrices F and F_opp.

Example files for 'mutation_data_breast.txt', 'breastDEgenes.txt' and 'network_genes.txt' are given.

'RANK_AGGREG' determines the rank aggregation scheme to be used. It can take values "SUM", "MEDIAN" or "RRA". "SUM" computes the summation of the per sample ranks and "MEDIAN" computes the median. "RRA" implements the robust rank aggregation technique as described in (Kolde et al, 2012). The matlab code for implementing the RRA method can be found at http://ch.mathworks.com/matlabcentral/fileexchange/41835-rank-aggregation. You will need to include the files betaScores.m, correctBetaPvalues.m, rhoScores.m and thresholdBetaScore.m.

'UP_DOWN' determines the diffusion procedure to be used. Possible entries are 1, 2 and 3. "1" performs diffusion from the genetically aberrant genes towards the directionality of the network interactions, "2" performs diffusion from the differentially expressed genes opposite from the directionality of the network interactions, and "3" performs bidirectional network diffusion (from both genetically aberrant and differentially expressed genes).

F and F_opp are precomputed diffused matrices on a given network. Given a directed adjacency matrix _adj_, the _F_ matrix can be computed as:

```
>> W = norm_adj( adj );
>> F = insulated_diff( W, .5 );
```

We first need to compute the normalized adjaceny matrix _W_ with the function _norm_adj_. The example above computes the _F_ matrix for a given restart probability 0.5. The adjacency matrix for the network described in (Wu et al, 2010) is given as a .mat file (_adj_lar_com.mat_). You can load the adjacency matrix by typing:

```
>> load('adj_lar_com');
```

The file _adj_lar_com_opp.mat_ is the transpose of the adjacency matrix. It can be used to compute _F_opp_ in the same way.

After the execution of the _netics_fun_, we can access the 10 highest ranked genes of the method by typing:

```
>> ranked_list_genes(1:10)
```

Dependencies:
  - Matlab (at least R2015a)

### Contributions
- [Christos Dimitrakopoulos](https://www.bsse.ethz.ch/cbg/group/people/person-detail.html?persid=197642)
- [Niko Beerenwinkel](http://www.bsse.ethz.ch/cbg/group/people/person-detail.html?persid=149417)


###Contact
```
Christos Dimitrakopoulos
christos.dimitrakopoulos (at) bsse.ethz.ch
```
