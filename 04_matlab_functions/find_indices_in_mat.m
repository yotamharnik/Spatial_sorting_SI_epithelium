function inds = find_indices_in_mat(mat_genes,genes)
% this fucntion will get a mat and will find the indices of the genes
% input:            mat: rows - gene exp values; coulmns - zones.
%                       exp values shoued be fraction of mRNA (devided by
%                       the sum)
%                   genes: list of genes that we want to find

inds = [];
for i = 1:length(genes)
    ind = find(strcmpi(mat_genes,genes(i)));
    if ~isempty(ind)
    inds = [inds ind(1)];
    end
end
