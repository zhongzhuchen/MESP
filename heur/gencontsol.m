function [contsol] = gencontsol(c,s)
n=size(c,1);
try
    [v,e]=eig(c);
catch
    [v,e]=eigs(c,n);
end
w=[diag(e),transpose(v)];
sw=sortrows(w,[1]);
v=transpose(sw(:,2:n+1));
contsol=zeros(n,1);

for j=1:n
  for i=1:s
    contsol(j)=contsol(j)+v(j,n-s+i)^2;
  end;
end;
