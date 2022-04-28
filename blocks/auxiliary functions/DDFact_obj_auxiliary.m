function [obj,dx,info] = DDFact_obj_auxiliary(x,s,F,Fsquare)
% This function calculate the objective value and gradient of DDFact for given data
%{
Input:
x       - current point for the DDFact problem
s       - the size of subset we want to choose, also equals to the summation of all elements of x
F       - C=FF' is a factorization of C where F is an n-by-d array
Fsquare - a 3d array where Fsquare(:,:,i) represents the F(i,:)'*F(i,:)

Output:
obj     - objective value of DDFact at current point x
dx      - the supgradient of obejctive function of DDFact at x
info    - recording of important information
        prob_nonsmooth - indicator of possible nonsmooth point
        dualgap - duality gap at current point x, which should be zero if x is an
                  optimal solution to DDFact
        dual_v  - corresponding dual solution
        dual_nu - corresponding dual solution
%}
n=length(x);
d=length(Fsquare(:,:,1));
info=struct;

X=zeros(d);
for i=1:n
    if x(i)==0
    else
        X=X+x(i)*Fsquare(:,:,i);
    end
end

% vectorized code of the above commented code, the experiments showed this
% vectorized version does not help

% xind=x>1e-12;
% x3d=zeros(1,1,n);
% x3d(1,1,:)=x;
% X=sum(Fsquare(:,:,xind).*x3d(:,:,xind),3);

X=1/2*(X+X');
[U,D]=eig(X);
D=diag(D);
sort_D=sort(D, 'descend');

% calculate k and corresponding mid_value as shown in Nikolov's paper
[k,mid_val]=find_k(sort_D,s);

if mid_val<=0
    sprintf('k=%d, mid_val=%f, sort_D(s)=%f, sum(x)=%f, rank(X)=%d', k, mid_val, sort_D(s), sum(x), rank(X))
    error('Something went wrong with calculating X or C might be a zero matrix.');
end

if abs(sort_D(k+1)-mid_val)<1e-6
    info.prob_nonsmooth=1;
else
    info.prob_nonsmooth=0;
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
ind1=D>mid_val;
eigDual(ind1)=1./D(ind1);
if mid_val>0
    eigDual(~ind1)=1/mid_val;
else
    eigDual(~ind1)=1/sort_D(k);
end


% calculate supgradient dx
% dX=U*diag(eigDual)*U';
% dx=diag(F*dX*F');

% vectorized code of the above commented code
K1=F*U;
K2=K1.*(eigDual');
dx=sum(K2.*K1,2);

% calculate objective value
sort_eigDual=sort(eigDual);
obj=-sum(log(sort_eigDual(1:s)));
end

