function [ind,val] = greedy(C,s)
% greedy algorithm without linear constraints
n=size(C,1);
[~,ind]=max(diag(C));
compind=setdiff(transpose(1:n),ind);
i=1;
while i<s
    i=i+1;
    j=0;
    candval=-10000000000;
    while j<length(compind)
        j=j+1;
        k=compind(j);
        h=log(det(C(union(ind,k),union(ind,k))));
        if h > candval
          cand=k;
          candval=h;
        end
    end
    ind=union(ind,cand);
    compind=setdiff(transpose(1:n),ind);
end
val=log(det(C(ind,ind)));
end