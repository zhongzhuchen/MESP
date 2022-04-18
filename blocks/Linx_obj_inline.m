%% obtain class properties and assign values
% scale C with Gamma
n = obj.size;
A_data=obj.A;
b_data=obj.b;
[m,~] = size(obj.A);
C=diag(Gamma)*obj.C;

%% calculate the objective value and gradient
F=C*diag(x)*C'+eye(n)-diag(x);
F=0.5*(F+F');

[R,flag]=chol(F); % F=R'*R cholesky decomposition

if flag>0
    warning("F(x) is not positive definite when calculating linx bound objective function.");
    fval=-Inf;
    dx=zeros(n,1);
    info.dual_upsilon=nan*ones(n,1);
    info.dual_nu=nan*ones(n,1);
    info.dualgap=nan;
    info.dualbound=nan;
else
    fval=sum(log(diag(R)))-sum(x.*log(Gamma)); % calculate the objective function
    Rinv=inv(R);
    K=C'*Rinv;
    % calculate the derivative: 1/2*diag(C'*F^{-1}*C-F^{-1})
    dx2=sum(Rinv.*Rinv,2);
    dx=0.5*(sum(K.*K,2)-dx2)-log(Gamma);
    
    %% build dual solution
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
    info.dualgap=dualgap+0.5*sum(dx2)-n/2+sum(x.*log(Gamma));
    info.fval=fval;
    info.dualbound=fval+info.dualgap;
    % cache for mixing
    info.cache = 0.5*sum(dx2)-n/2;
end


