# NetICS
Network-based integration of multi omics data for prioritizing cancer genes

NetICS performs a per sample bidirectional network diffusion technique to prioritize genes based on their proximity to genetically aberrant and differentially expressed genes. It provides a rank aggregation technique for integrating the sample-specific gene lists into an overall ranked list of genes.

The method should be called as follows:

```
ranked_list_genes = netics_fun( breast_struct, F_04, F_04_opp, 'RANK_AGGREG=SUM', 'UP_DOWN=3', 'breastDEgenes.txt', 'network_genes.txt');

```

'breastDEgenes.txt' --> a file containing a list of gene names detected as differentially expressed.
'network_genes.txt' --> a file containing the list of network gene names. They should be in the same order as in the diffused matrices F_04 and F_04_opp.

Example files for 'breastDEgenes.txt' and 'network_genes.txt' are given.

RANK_AGGREG can take values "SUM" or "MEDIAN". "SUM" computes the summation of the per sample ranks and "MEDIAN" computes the median.

UP_DOWN can take values from 1 to 3. "1" performs diffusion from the genetically aberrant genes, "2" performs diffusion from the differentially expressed genes, and "3" performs bidirectional network diffusion (from both aberrant and differentially expressed genes).

### Contributions
- [Christos Dimitrakopoulos](https://www.bsse.ethz.ch/cbg/group/people/person-detail.html?persid=197642)
- [Niko Beerenwinkel](http://www.bsse.ethz.ch/cbg/group/people/person-detail.html?persid=149417)


###Contact
```
Christos Dimitrakopoulos
christos.dimitrakopoulos (at) bsse.ethz.ch
```
