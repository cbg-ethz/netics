function [ ranked_list_genes, scores ] = netics_fun( varargin )

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  
%   INPUTS:
%   filenameMu(varargin{1}):         input file that contains the genetically aberrant genes of each sample. 
%                                    It contains two columns that map every gene (1st column) to the samples that 
%                                    it is genetically aberrant (2nd column)
%   adj(varargin{2}):                adjacency matrix of the directed interaction network
%   restart_prob(varargin{3}):       restart probability for the insulated diffusion. For (Wu et al., 2010) network use 0.4.
%   rank_method_str(varargin{4}):    'MEDIAN' uses the median of the sample-specific ranks
%                                    'RRA' uses the Robust Rank Aggregation method to integrate sample-specific ranked lists
%                                    'SUM' uses the sum of the sample-specific ranks
%   filenameNet(varargin{5}):        input file that contains the list of the genes that are present in the network. 
%                                    They should be in the same order as in the rows of the adjacency matrix adj. 
%                                    An example file is given that contains the gene names of the network described in (Wu et al, 2010).
%   filenameRNA(varargin{6}):        tab delimited file with two columns. First column contains the genes for which differential expression 
%                                    between the tumor and normal samples at the RNA level was measured. Second column contains the p-values 
%                                    of these measurements. This file can be the result of a tool for differential expression analysis such as DESeq2. 
%                                    Each gene in this file should have only one entry.
%   filenamePR(varargin{7}):         tab delimited file with two columns. First column contains the proteins for which differential expression between 
%                                    the tumor and normal samples at the protein level was measured. Second column contains the p-values of these 
%                                    measurements. Each gene in this file should have only one entry.
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  
%   OUTPUTS:
%   ranked_list_genes:  list with the ranked gene names
%   scores:             list with the scores of the ranked gene names
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  

	filenameRNA = '';
	filenamePR = '';
	ranked_list_genes={};

	% Accept arguments    
	if nargin < 5,
		error('Error: Missing input arguments.')
	elseif nargin > 7,
		error('Error: Too many input arguments.')
	else            
		mutation_data = read_mutations( varargin{1} );
		adj = varargin{2};
		restart_prob = varargin{3};
		rank_method_str = varargin{4};
		filenameNet = varargin{5};
		choose_mut_diff = 1;
	end
	if nargin >= 6,            
		filenameRNA = varargin{6};
		choose_mut_diff = 3;        
	end
	if nargin >= 7,
		filenamePR = varargin{7};
	end

	%read network genes
	fidnet = fopen(filenameNet,'r');
	g = textscan(fidnet, '%s', 'delimiter', '\n');
	network_genes = g{1};
	fclose(fidnet);
	
	%Extract first token
	rank_method_str1 = rank_method_str(1:strfind(rank_method_str,'='));
	
	%Handle input
	rank_method_str = rank_method_str(strfind(rank_method_str,'=')+1:end);
	
	%Check input - Token 1
	if ~strcmp(rank_method_str1,'RANK_AGGREG='),
		disp(strcat('Wrong input: ',rank_method_str1));
		disp('Input should be RANK_AGGREG=');		
		return;	
	end

	%Check input - Token 2
	if  ~strcmp(rank_method_str,'SUM') && ~strcmp(rank_method_str,'MEDIAN') && ~strcmp(rank_method_str,'RRA'),
		disp(strcat('Wrong input: ',rank_method_str));
		disp(strcat('Input should be MEDIAN, RRA or SUM'));
		return;		
	end

	%copy data from input mutation_data 
	samples = {};
	for i = 1:length(mutation_data),
		samples.DNA{i} = upper(mutation_data{i});
		samples.DNA_network{i} = upper(intersect(network_genes,samples.DNA{i}));
	end
	clear mutation_data;

	%read differentially expressed genes
	diff_expr_genes = {};
	if choose_mut_diff~=1,
		diff_expr_genes = read_diff_expr( filenameRNA, filenamePR );
		if isempty(diff_expr_genes),
			choose_mut_diff=1;
		end
	else
		disp('No Input files were provided for differentially expressed genes.');        
	end
	for i = 1:length(samples.DNA),
		samples.RNA{i} = diff_expr_genes;
		samples.RNA_network{i} = upper(intersect(network_genes,samples.RNA{i}));
	end

	%compute diffused matrix
	disp('Computing diffused matrix...');
	F = insulated_diff(norm_adj(adj), restart_prob);
	F_opp = insulated_diff(norm_adj(adj'), restart_prob);    
	disp('Running NetICS...');
	[ranked_list_genes, scores] = prioritization( samples, F, F_opp, network_genes, choose_mut_diff, rank_method_str);

end

function [ mutation_data ] = read_mutations( filename )
	fid = fopen(filename,'r');
	g = textscan(fid, '%s%s', 'delimiter', '\t');
	fclose(fid);
	unique_samples = unique(g{2});
	mutation_data = cell(length(unique_samples),1);
	for i = 1:length(unique_samples),
		mutation_data{i} = g{1}(ismember(g{2},unique_samples{i}));
	end
end

function [ diff_expr_genes ] = read_diff_expr( filenameDE, filenamePR )

	diff_expr_genes ={}; RNA_names={}; rppa_names={};

	if ~isempty(filenamePR),
		[rppa_names, rppa_pval] = read_file(filenamePR);
	else
		disp('No Input file was provided for differentially expressed genes at the proteome level.');        
	end

	if ~isempty(filenameDE),
		[RNA_names, RNA_pval] = read_file(filenameDE);        
	else
		disp('No Input file was provided for differentially expressed genes at the RNA level.');        
	end

	%check for duplicate genes in the input files
	if (length(unique(RNA_names)) ~= length(RNA_names)) || (length(unique(rppa_names)) ~= length(rppa_names)),
		disp('Input files for differentially expressed genes should contain only one entry per gene.');
		disp('Input files for differentially expressed genes were ignored.');        
		return;	        
	end
	%check if the input files for differential expression are empty    
	if isempty(rppa_names) && isempty(RNA_names),
		disp('Only genetically aberrant genes are used for diffusion.');
		return;
	%only differential expressed genes at the RNA level are provided
	elseif isempty(rppa_names) && ~isempty(RNA_names),
		[~,~,pval_all] = fdr(RNA_pval,0.05);
		diff_expr_genes = RNA_names(pval_all < 0.05);
	%only differential expressed genes at the protein level are provided        
	elseif isempty(RNA_names) && ~isempty(rppa_names),
		[~,~,pval_all] = fdr(rppa_pval,0.05);
		diff_expr_genes = rppa_names(pval_all < 0.05);
	else
		names_intersection = intersect(RNA_names,rppa_names);
		names_setdiff_rna = setdiff(RNA_names,rppa_names);
		names_setdiff_pr = setdiff(rppa_names,RNA_names);

		rna_pval_names_intersection = RNA_pval(ismember(RNA_names,names_intersection));
		rppa_pval_names_intersection = rppa_pval(ismember(rppa_names,names_intersection));

		rna_pval_names_setdiff_rna = RNA_pval(ismember(RNA_names,names_setdiff_rna));
		pr_pval_names_setdiff_pr = rppa_pval(ismember(rppa_names,names_setdiff_pr));

		all_names = [names_intersection' names_setdiff_rna' names_setdiff_pr'];
		%Fishers method
		pval_all = (1-pchisq(-2*(log(rna_pval_names_intersection)+log(rppa_pval_names_intersection)),4));
		pval_all = [pval_all' rna_pval_names_setdiff_rna' pr_pval_names_setdiff_pr'];    

		[~,~,pval_all] = fdr(pval_all,0.05);
		diff_expr_genes = all_names(pval_all < 0.05);
	end
end

function [names, pval] = read_file(filenamePR)
	fid = fopen(filenamePR,'r');
	g = textscan(fid, '%s%s', 'delimiter', '\t');
	fclose(fid);
	names = g{1}; pval = str2double(g{2}); 
end
