function high_protein_ind = find_most_high_exp_protein_in_mrna(proteins,st)
% this function gets a list of proteins names and a structure of RNA data.
% it looks for what is the most highly expressed protein in the mrna level
% and gives his index in the protein list
%   ** if each of the proteins is not found in the mRNA struct, the
%   function will return the 1st index as the high protein
% input:     proteins        >> a list of ptotein names
%            st.mat_norm     >> genes over samples mat (usally norm by the sum)
%            st.mean_norm    >> genes over zones mat (mean over all mice)
%            st.gene_name    >> gene names in the same order of the roes in the mat
%            st.sample_names >> holds the sample name in the same order in
%                               the mat for the plotting           
%  output:   high_protein_ind >> index of the most highly expressed protein

NO_MATCH_FLAG = 1;
exp_mat = zeros(length(proteins),1);
for i = 1 : length(proteins)
    ind = find(strcmpi(st.gene_name,proteins(i)));
    if ~isempty(ind)
        exp_mat(i) = mean(st.mean_norm(ind(1),:));
        NO_MATCH_FLAG = 0;
    end
end
if NO_MATCH_FLAG
    high_protein_ind = 1;
else
[~,high_protein_ind] = max(exp_mat);
end