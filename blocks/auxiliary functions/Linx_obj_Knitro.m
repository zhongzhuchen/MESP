function [fval,dx] = Linx_obj_Knitro(x,C,Gamma)
% light version of Linx_obj without 
%% obtain class properties and assign values
% scale C with Gamma
n = length(x);
C=diag(Gamma)*C;

%% calculate the objective value and gradient
F=C*diag(x)*C'+eye(n)-diag(x);
F=0.5*(F+F');

[R,flag]=chol(F); % F=R'*R cholesky decomposition

if flag>0
    warning("F(x) is not positive definite when calculating linx bound objective function.");
    fval=-Inf;
    dx=zeros(n,1);
else
    fval=sum(log(diag(R)))-sum(x.*log(Gamma)); % calculate the objective function
    Rinv=inv(R);
    K=C'*Rinv;
    % calculate the derivative: 1/2*diag(C'*F^{-1}*C-F^{-1})
    dx2=sum(Rinv.*Rinv,2);
    dx=0.5*(sum(K.*K,2)-dx2)-log(Gamma);
end

fval=-fval;
dx=-dx;
end

