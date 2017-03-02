# NetICS: Network-based integration of multi-omics data for prioritizing cancer genes

NetICS performs a per sample bidirectional network diffusion technique to prioritize genes based on their proximity to genetically aberrant and differentially expressed genes. It provides rank aggregation techniques for integrating the sample-specific gene lists into an overall ranked list of genes.

The method is called as follows:

```
>> [ ranked_list_genes, ranked_scores ] = netics_fun( 'mutation_data_breast.txt', F, F_opp, 'RANK_AGGREG=SUM', 'network_genes.txt', 'RNA_diff_expr_breast.txt', 'protein_diff_expr_breast.txt');
```

'mutation_data_breast.txt' --> tab delimited file that contains the genetically aberrant genes of each sample. It contains two columns that map each gene (1<sup>st</sup> column) to the samples that it is genetically aberrant (2<sup>nd</sup> column).

'network_genes.txt' --> input file that contains the list of the genes that are present in the network. They should be in the same order as in the rows of the diffused matrices F and F_opp. An example file is given that contains the gene names of the network described in (Wu et al, 2010).

'RNA_diff_expr_breast.txt' --> tab delimited file with two columns. First column contains the genes for which differential expression between the tumor and normal samples at the RNA level was measured. Second column contains the p-values of these measurements. This file can be the result of a tool for differential expression analysis such as DESeq2. Each gene should have only one entry in this file.

'protein_diff_expr_breast.txt' --> tab delimited file with two columns. First column contains the proteins for which differential expression between the tumor and normal samples at the protein level was measured. Second column contains the p-values of these measurements. Each gene should have only one entry in this file.

The two files that contain the differentially expressed genes at the RNA and proteome levels (for example, 'RNA_diff_expr_breast.txt' and 'protein_diff_expr_breast.txt') are optional. If not provided, only the labels of the genetically aberrant genes are used for network diffusion.

The p-values in files 'RNA_diff_expr_breast.txt' and 'protein_diff_expr_breast.txt' should be provided unadjusted because they are combined by using the Fisher's method. After that, NetICS adjusts them for multiple testing by using Benjamini & Hochberg FDR correction. The function for computing FDR correction can be derived from https://brainder.org/2011/09/05/fdr-corrected-fdr-adjusted-p-values/.

Example files for 'mutation_data_breast.txt', 'RNA_diff_expr_breast.txt' and protein_diff_expr_breast.txt' are given based on the breast invasive carcinoma dataset in TCGA (BRCA).

'RANK_AGGREG' determines the rank aggregation scheme to be used. It can take the values "SUM", "MEDIAN" or "RRA". "SUM" computes the summation of the per sample ranks and "MEDIAN" computes the median. "RRA" implements the robust rank aggregation technique as described in (Kolde et al, 2012). The matlab code of the RRA method can be derived from http://ch.mathworks.com/matlabcentral/fileexchange/41835-rank-aggregation. You will need to include the files betaScores.m, correctBetaPvalues.m, rhoScores.m and thresholdBetaScore.m.

F and F_opp are precomputed diffused matrices on a given network. Given a directed adjacency matrix _adj_, the _F_ matrix can be computed as:

```
>> W = norm_adj( adj );
>> F = insulated_diff( W, .5 );
```

We first need to compute the normalized adjaceny matrix _W_ with the function _norm_adj_. The above example computes the _F_ matrix for a given restart probability 0.5. The adjacency matrix for the network described in (Wu et al, 2010) is given as a .mat file (_adj_lar_com.mat_). You can load the adjacency matrix by typing:

```
>> load('adj_lar_com');
```

The file _adj_lar_com_opp.mat_ is the transpose of the adjacency matrix. It can be used to compute _F_opp_ in the same way.

After the execution of the _netics_fun_, we can access the 10 highest ranked genes of the method by typing:

```
>> ranked_list_genes(1:10)
```
The files _pchisq.m_ and _pgamma.m_ were derived from https://ch.mathworks.com/matlabcentral/fileexchange/15171-jennrich-test/content/Jennrich/pchisq.m.

Whenever the word 'sample' is mentioned above, we mean one paired observation for which measurements for tumor and normal tissues are available.

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
