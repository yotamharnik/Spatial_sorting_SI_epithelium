function [com_vec] = calculate_com_mat(mat,com_lim)
if nargin == 1
    com_lim = 1;
end
com_vec = zeros(size(mat,1),1);
for i = 1 :size(mat,1)
    vec=(mat(i,:));
    com_vec(i) = nansum(linspace(0,com_lim,length(vec)).*vec)/nansum(vec);
end
end