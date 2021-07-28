function A = read_gene_list(file)
% this function get a txt file and produces a gene names list.
% input:    path to a txt file. one gene name in a line
% output:   a cell array of gene names without duplicates/
fileID = fopen(file);
B = textscan(fileID,'%s');
A = unique(B{1});
end
