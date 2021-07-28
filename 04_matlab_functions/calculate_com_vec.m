function [com] = calculate_com_vec(vec,com_lim)
if nargin == 1
    com_lim = 1;
end
   com = nansum(linspace(0,com_lim,length(vec)).*vec)/nansum(vec); 
end