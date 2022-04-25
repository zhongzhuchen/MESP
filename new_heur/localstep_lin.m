function [xind,val] = localstep_lin(C,xind,A,b)
% local search under linear constraints
n=size(C,1);
s=length(xind);
compind=setdiff(transpose(1:n),xind);
val=log(det(C(xind,xind)));
i=0;
while i<s
  i=i+1;
  j=0;
  while j < n-s
    j=j+1;
    k=xind(i);
    l=compind(j);
    newxind = union(setdiff(xind,k),l);
    newval=log(det(C(newxind,newxind)));
    newx = zeros(n,1);
    newx(newxind)=1;

    if newval > val && sum(max(A*newx-b,0))<=1e-12
        xind=newxind;
        compind=union(setdiff(compind,l),k);
        val=newval;
        i=0;
        break;
    end

  end
end
end

