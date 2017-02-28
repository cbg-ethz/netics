function [ ranked_list_genes, scores ] = netics_fun( filenameMu, F, F_opp, rank_method_str, choose_mut_diff_flag, filenameRNA, filenamePR, filenameNet)

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  
%   INPUTS:
%   filenameMu:         input file that contains the genetically aberrant genes of each sample. 
%                       It contains two columns that map every gene (1st column) to the samples that 
%                       it is genetically aberrant (2nd column)
%   F:                  diffused network towards the directionality of the network interactions
%   F_opp:              diffused network opposite of the directionality of the network interactions
%   rank_method_str:    'MEDIAN' uses the median of the sample-specific ranks
%                       'RRA' uses the Robust Rank Aggregation method to integrate sample-specific ranked lists
%                       'SUM' uses the sum of the sample-specific ranks
%   choose_up_down_flag: ==1 diffuse only mutations, ==2 diffuse only RNA, ==3 diffuse both
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  
%   OUTPUTS:
%   ranked_list_genes:  list with the ranked gene names
%   scores:             list with the scores of the ranked gene names
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  

	ranked_list_genes={};
	
	%read mutation data
	mutation_data = read_mutations( filenameMu );

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
	diff_expr_genes = {};
	if choose_mut_diff~=1,
		diff_expr_genes = read_diff_expr( filenameRNA, filenamePR );
		if isempty(diff_expr_genes),
			return;
		end
	end
	for i = 1:length(samples.DNA),
		samples.RNA{i} = diff_expr_genes;
		samples.RNA_network{i} = upper(intersect(network_genes,samples.RNA{i}));
	end 
	[ranked_list_genes, scores] = prioritization( samples, F, F_opp, network_genes, choose_mut_diff, rank_method_str);

end

function [ mutation_data ] = read_mutations( filename )
	fid = fopen(filename,'r');
	g = textscan(fid, '%s%s', 'delimiter', '\t');
	fclose(fid);
	unique_samples = unique(g{2});
	for i = 1:length(unique_samples),
		pos = ismember(g{2},unique_samples{i});
		mutation_data{i}.genes = g{1}(pos);
	end
end

function [ diff_expr_genes ] = read_diff_expr( filenameRNA, filenamePR )
	diff_expr_genes ={};    
	fid = fopen(filenameRNA,'r');
	gRNA = textscan(fid, '%s%s', 'delimiter', '\t');
	fclose(fid);
	fid = fopen(filenamePR,'r');
	gPR = textscan(fid, '%s%s', 'delimiter', '\t');
	fclose(fid);
	rppa_names = gPR{1}; rppa_pval = str2double(gPR{2});
	RNA_names = gRNA{1}; RNA_pval = str2double(gRNA{2});
	clear gRNA gPR;

	%check for duplicate genes in the input files
	if (length(unique(RNA_names)) ~= length(RNA_names)) || (length(unique(rppa_names)) ~= length(rppa_names)),
		disp('Input files for differentially expressed genes should contain only one entry per gene.');
		disp('Input files for differentially expressed genes were ignored.');        
		return;	        
	end
	%check if the input files for differential expression are empty    
	if isempty(rppa_names) && isempty(RNA_names),
		disp('Input files for differentially expressed genes should be provided!');
		return;
	%only differential expressed genes at the RNA level are provided
	elseif isempty(rppa_names) && ~isempty(RNA_names),
		[~,~,pval_all] = fdr(RNA_pval,0.05);
		diff_expr_genes = RNA_names(pval_all < 0.05);
		disp('Differential expression at the protein level was not provided.');
	%only differential expressed genes at the protein level are provided        
	elseif isempty(RNA_names) && ~isempty(rppa_names),
		[~,~,pval_all] = fdr(rppa_pval,0.05);
		diff_expr_genes = rppa_names(pval_all < 0.05);
		disp('Differential expression at the RNA level was not provided.');
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
