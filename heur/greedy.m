function [ind,val] = greedy(c,n,s)
[a,b]=max(diag(c));
ind=b;
compind=setdiff(transpose(1:n),ind);
i=1;
while i<s,
  i=i+1;
  j=0;
  candval=-10000000000;
  while j<n-s,
    j=j+1;
    k=compind(j);
    h=log(det(c(union(ind,k),union(ind,k))));
    if h > candval
      cand=k;
      candval=h;
    end;
  end;
ind=union(ind,cand);
compind=setdiff(transpose(1:n),ind);
end;
val=log(det(c(ind,ind)));
