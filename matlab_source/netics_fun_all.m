%script that prioritize genes by taking the mean ranked position across several values of the restart probability
%store final ranked list of genes in ranked_list_genes

fidnet = fopen('network_genes.txt','r');
g = textscan(fidnet, '%s', 'delimiter', '\n');
network_genes = g{1};
fclose(fidnet);

load('adj_lar_com');
betas = .2:.1:.8;

ranked_list_genes_breast = cell(length(betas),1);
for i = 1:length(betas),
    disp(strcat('Restart Probability = ',num2str(betas(i))));
    ranked_list_genes_breast{i} = netics_fun( 'mutation_data_breast.txt', adj_lar_com, betas(i), 'RANK_AGGREG=SUM', 'network_genes.txt', 'RNA_diff_expr_breast.txt', 'protein_diff_expr_breast.txt');
end

ranks = zeros(1,length(network_genes));
for i = 1:length(network_genes),
    pos = zeros(1,length(betas));
    for j = 1:length(betas),
        pos(j) = find(ismember(ranked_list_genes_breast{j},network_genes{i}));
    end
    ranks(i) = mean(pos);
end

sorted = sortrows([ranks; 1:length(ranks)]');
ranked_list_genes = network_genes(sorted(:,2));
