function [contsol] = gencontsol_eig(C,s)
n=size(C,1);
try
    [v,e]=eig(C);
catch
    [v,e]=eigs(C,n);
end
w=[diag(e),transpose(v)];
sw=sortrows(w,[1]);
v=transpose(sw(:,2:n+1));

contsol = sum(v(:, (n-s+1):n).^2, 2);
% contsol1=zeros(n,1);
% 
% for j=1:n
%   for i=1:s
%     contsol1(j)=contsol1(j)+v(j,n-s+i)^2;
%   end
% end
end
