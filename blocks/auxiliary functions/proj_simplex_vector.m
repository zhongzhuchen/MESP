function [y] = proj_simplex_vector(y)
% credit to Laurent Condat @https://lcondat.github.io/software.html
y=max(y-max((cumsum(sort(y,1,'descend'),1)-1)./(1:size(y,1))'),0);
end

