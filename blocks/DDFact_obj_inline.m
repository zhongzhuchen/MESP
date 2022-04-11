%% obtain class properties and assign values
% scale F, Fsquare with Gamma
n=obj.size;
F=diag(sqrt(Gamma))*obj.F;
Fsquare = obj.Fsquare;
for i=1:n
    Fsquare(:,:,i)=Gamma(i)*Fsquare(:,:,i);
end
n = obj.size;
A_data=obj.A;
b_data=obj.b;
[m,~] = size(obj.A);
%% obtain parameter properties
d=length(Fsquare(:,:,1));

%% a struct containing information
info=struct;

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
%     obj=-Inf;
%     dx=zeros(n,1);
%     info=0;
%     sprintf('k=%d, mid_val=%f, sort_D(s)=%f, rank(X)=%d', k, mid_val, sort_D(s), rank(X))
%     return
    sprintf('k=%d, mid_val=%f, sort_D(s)=%f, sum(x)=%f, rank(X)=%d', k, mid_val, sort_D(s), sum(x), rank(X))
    error('Something went wrong with calculating X or C might be a zero matrix.');
end

% construct the eigenvalue for the feasible solution to the dual problem, i.e., DFact
eigDual=zeros(d,1);
% for i=1:d
%     if(D(i)>mid_val)
%         eigDual(i)=1/D(i);
%     else
%         if mid_val>0
%             eigDual(i)=1/mid_val;
%         else
%             eigDual(i)=1/sort_D(k);
%         end
%     end
% end

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
%% calculate dual bound and dual solutions by calling Knitro
% % cache for mixing
% info.cache1=-s;
% info.cache2=sum(dx)-s;

f=[zeros(n,1);ones(n,1);b_data;s];
Aeq=[-eye(n),eye(n),A_data',ones(n,1)];
beq=dx;
lb=[zeros(2*n+m,1);-inf];
ub=Inf(2*n+m+1,1);
x0=[];
options = knitro_options('algorithm',3,...  % active-set/simplex algorithm
                         'outlev',0);       % iteration display
[xlp, dualgap, exitflag, ~] = knitro_lp (f, [], [], Aeq, beq, lb, ub, x0, [], options);
info.dual_upsilon = xlp(1:n);
info.dual_nu = xlp((n+1):2*n);
info.dual_pi =  xlp((2*n+1):(2*n+m));
info.dual_tau = xlp(end);
info.dualgap=dualgap-s+sum(x.*log(Gamma));
% fval=-sum(log(eigDual(ind1)));
sort_eigDual=sort(eigDual);
fval=-sum(log(sort_eigDual(1:s)))-sum(x.*log(Gamma));
info.fval=fval;
info.dualbound=fval+info.dualgap;
