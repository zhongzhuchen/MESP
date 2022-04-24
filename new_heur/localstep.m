function [xind,val] = localstep(C,xind)
% local search without linear constraints
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
    if newval > val
        xind=newxind;
        compind=union(setdiff(compind,l),k);
        val=newval;
        i=0;
        break;
    end
  end
end
end
