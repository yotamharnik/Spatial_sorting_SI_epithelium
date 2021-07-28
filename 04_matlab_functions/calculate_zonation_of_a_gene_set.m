function [zonation_vec,zonation_sem] = calculate_zonation_of_a_gene_set(mat)
% this function will get a mat of n genes over m zones and will preduce an
% avraged zonation profile.
% first all the genes will be normlised to the mean value (so they will be in
% same sacle) and then we will take from each zone the median expression.
mat_norm = mat./nanmean(mat,2);
zonation_vec = nanmean(mat_norm);
zonation_sem = nanstd(mat_norm)/sqrt(size(mat_norm,1));
end