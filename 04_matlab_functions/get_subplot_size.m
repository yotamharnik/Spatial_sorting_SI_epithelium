function [m,n] = get_subplot_size(gg)
    m   =   ceil(sqrt(length(gg)));
    n   =   ceil(length(gg)/m);
end