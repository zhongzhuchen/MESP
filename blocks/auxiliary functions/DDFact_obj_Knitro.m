function [fval,dx] = DDFact_obj_Knitro(x,s,F,Fsquare,Gamma)
%% obtain class properties and assign values
% scale F, Fsquare with Gamma
n=length(x);
F=diag(sqrt(Gamma))*F;
for i=1:n
    Fsquare(:,:,i)=Gamma(i)*Fsquare(:,:,i);
end
%% obtain parameter properties
d=length(Fsquare(:,:,1));

X=zeros(d);
for i=1:n
    if x(i)==0
    else
        X=X+x(i)*Fsquare(:,:,i);
    end
end

X=1/2*(X+X');
[U,D]=eig(X);
D=diag(D);
[sort_D, sort_idx]=sort(D, 'descend');

%% calculate k and corresponding mid_value as shown in Nikolov's paper
[k,mid_val]=find_k(sort_D,s);

if mid_val<=0
    sprintf('k=%d, mid_val=%f, sort_D(s)=%f, sum(x)=%f, rank(X)=%d', k, mid_val, sort_D(s), sum(x), rank(X))
    error('Something went wrong with calculating X or C might be a zero matrix.');
end

% construct the eigenvalue for the feasible solution to the dual problem, i.e., DFact
eigDual=zeros(d,1);

% vectorized code of the above commented code
ind1=sort_idx(1:k);
non_ind1= setdiff(1:d, ind1);
eigDual(ind1)=1./D(ind1);
if mid_val>0
    eigDual(non_ind1)=1/mid_val;
else
    eigDual(non_ind1)=1/sort_D(k);
end

%% calculate the gradient dx
% dX=U*diag(eigDual)*U';
% dx=diag(F*dX*F');

% vectorized code of the above commented code
K1=F*U;
K2=K1.*(eigDual');
dx1=sum(K2.*K1,2);
dx=dx1-log(Gamma);
sort_eigDual=sort(eigDual);
fval=-sum(log(sort_eigDual(1:s)))-sum(x.*log(Gamma));

fval=-fval;
dx=-dx;
end
