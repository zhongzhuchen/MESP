function [xind,val] = localstep(c,xind)
n=size(c,1);
s=size(xind);
compind=setdiff(transpose(1:n),xind);
val=log(det(c(xind,xind)));
i=0;
j=0;
while i<s,
  i=i+1;
  while j < n-s,
    j=j+1;
    k=xind(i);
    l=compind(j);
    valnew=log(det(c(union(setdiff(xind,k),l),union(setdiff(xind,k),l))));
    if valnew > val
    xind=union(setdiff(xind,k),l);
    compind=union(setdiff(compind,l),k);
    val=valnew;
    i=0;
    j=0;
    break;
  end;
end;

end;
