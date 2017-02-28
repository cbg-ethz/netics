function [ ranked_list_genes, ranked_scores ] = prioritization( samples, F, F_opp, network_genes, choose_up_down_flag, rank_method_str)

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
%   INPUTS:
%   samples:                cell array of samples, contains the DNA changes and differentially expressed genes of each sample
%   F:                      diffused network (insulated heat diffusion).
%   F_opp:                  diffused network opposite direction (insulated heat diffusion).
%   network_genes:          the names of all the network genes in the same order they appear in the adjacency matrix.
%   choose_up_down_flag:    1 diffuse only from mutations, 2 diffuse only from differentially expressed genes, 3 diffuse from both.
%   rank_method_str:        'MEDIAN' uses the median of the sample-specific ranks.
%                           'RRA' uses the Robust Rank Aggregation method to integrate sample-specific ranked lists.
%                           'SUM' uses the summation of the sample-specific ranks.
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
%   OUTPUTS:
%   ranked_list_genes:      list with the ranked gene names
%   ranked_scores:          list with the scores of the ranked gene names
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
	
	weights_all = zeros(length(samples.DNA_network),length(network_genes));
	positions_sorted = zeros(length(network_genes),length(samples.DNA_network));
	[final_weights_median,final_weights_sum,final_weights_rho] = deal(zeros(1,length(network_genes)));

	for i = 1:length(samples.DNA_network), %number of samples
		weights_all(i,:) = diffuse_all(choose_up_down_flag, samples.DNA_network{i}, samples.RNA_network{i}, network_genes, F, F_opp);
		to_sort_sorted = sortrows([weights_all(i,:) ;1:length(weights_all(i,:))]');
		positions_sorted(:,i) = fliplr(to_sort_sorted(:,2)')';		
	end
	for j = 1:length(network_genes), % number of genes
		rank = zeros(length(samples.DNA_network),1);
		for i = 1:length(samples.DNA_network), % number of samples
			rank(i) = find(positions_sorted(:,i) == j);
		end
		final_weights_median(j) = median(rank);
		final_weights_sum(j) = sum(rank);
		final_weights_rho(j) = rhoScores(rank/max(rank));    
	end
	
	%choose the aggregation method 
	if strcmp(rank_method_str, 'MEDIAN'),
		scores = final_weights_median;
	elseif strcmp(rank_method_str, 'SUM'),
		scores = final_weights_sum;		
	else
		scores = final_weights_rho;
	end
	
	%sort scores
	to_sort_sorted_rho = sortrows([scores ;1:length(scores)]');
	ranked_list_genes = network_genes(to_sort_sorted_rho(:,2));
	ranked_scores = scores(to_sort_sorted_rho(:,2));
end

function [ret_weights_all] = diffuse_all(choose_up_down_flag, sample_DNA, sample_RNA, network_genes, F, F_opp)
	%diffuse mutation labels
	if choose_up_down_flag == 1 || choose_up_down_flag == 3,
		mutation_weights = diffuse(intersect(sample_DNA,network_genes), network_genes, F);
	end
	%diffuse DE labels	
	if choose_up_down_flag == 2 || choose_up_down_flag == 3,
		DE_weights = diffuse(intersect(sample_RNA,network_genes), network_genes, F_opp);
	end
	%combines scores by multiplication (use both mutations and DE)
	if choose_up_down_flag == 3,
		ret_weights_all = mutation_weights.*DE_weights;
	elseif choose_up_down_flag == 1,
		ret_weights_all = mutation_weights;
	else
		ret_weights_all = DE_weights;
	end
end	
	
function [mutation_weights] = diffuse(inter_mut_sampl_network, network_genes, F)
	positions = find(ismember(network_genes,inter_mut_sampl_network));
	mutation_weights= (1/length(positions))*ones(1,length(positions))*F(positions,:);
end
