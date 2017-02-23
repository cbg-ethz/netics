function [ ranked_list_genes ] = netics_fun( mutation_data, F, F_opp, rank_method_str, choose_mut_diff_flag, filenameDE, filenameNet)

	%INPUTS
	%mutation_data: cell array. Each cell contains a struct of data about the mutations/copy number of the samples. 
	%               access field "genes" for getting the aberrant genes of each sample.
	%F: diffused network (insulated heat diffusion).
	%F_opp: diffused network opposite direction (insulated heat diffusion).
	%rank_method_str: 'MEDIAN' uses the median of the sample-specific ranks.
	%				  'RRA' uses the Robust Rank Aggregation method to integrate sample-specific ranked lists.
	%				  'SUM' uses the sum of the sample-specific ranks.   
	%choose_up_down_flag: ==1 diffuse only mutations, ==2 diffuse only RNA, ==3 diffuse both
	%OUTPUTS
	%ranked_list_genes: list with the ranked gene names
	
	ranked_list_genes={};
	
	%read network genes
	fidnet = fopen(filenameNet,'r');
	g = textscan(fidnet, '%s', 'delimiter', '\n');
	network_genes = g{1};
	fclose(fidnet);
	
	%Extract first token
	rank_method_str1 = rank_method_str(1:strfind(rank_method_str,'='));
	choose_mut_diff1 = choose_mut_diff_flag(1:strfind(choose_mut_diff_flag,'='));
	
	%Handle input
	rank_method_str = rank_method_str(strfind(rank_method_str,'=')+1:end);
	choose_mut_diff = str2double(choose_mut_diff_flag(strfind(choose_mut_diff_flag,'=')+1:end));
	
	%Check input - Token 1
	if ~strcmp(rank_method_str1,'RANK_AGGREG='),
		disp(strcat('Wrong input: ',rank_method_str1));
		disp('Input should be RANK_AGGREG=');		
		return;	
	end
	if ~strcmp(choose_mut_diff1,'UP_DOWN='),
		disp(strcat('Wrong input: ',choose_mut_diff1));
		disp('Input should be UP_DOWN=');		
		return;	
	end	

	%Check input - Token 2
	if  ~strcmp(rank_method_str,'SUM') && ~strcmp(rank_method_str,'MEDIAN') && ~strcmp(rank_method_str,'RRA'),
		disp(strcat('Wrong input: ',rank_method_str));
		disp(strcat('Input should be MEDIAN, RRA or SUM'));
		return;		
	end
	if isempty(choose_mut_diff) || sum(choose_mut_diff~=[1 2 3])==3,
		disp(strcat('Wrong input: ',choose_mut_diff));
		disp(strcat('Input should be 1, 2, 3'));
		return;		
	end
    
	%copy data from input mutation_data 
	samples = {};
	for i = 1:length(mutation_data),
		samples.DNA{i} = upper(mutation_data{i}.genes);
		samples.DNA_network{i} = upper(intersect(network_genes,samples.DNA{i}));
	end
	
	%read differentially expressed genes
	fiddeseq = fopen(filenameDE,'r');
	g = textscan(fiddeseq, '%s', 'delimiter', '\n');
	RNA = unique(g{1});
	fclose(fiddeseq);
	for i = 1:length(samples.DNA),
		samples.RNA{i} = RNA;
		samples.RNA_network{i} = upper(intersect(network_genes,samples.RNA{i}));
	end
	ranked_list_genes = prioritization( samples, F, F_opp, network_genes, choose_mut_diff, rank_method_str);
	
end
