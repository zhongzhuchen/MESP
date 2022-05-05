%% obtain class properties
C=obj.C;
n = obj.size;
A_data=obj.A;
b_data=obj.b;
[m,~] = size(obj.A);
F=obj.F;
Fsquare=obj.Fsquare;
F_comp=obj.F_comp;
Fsquare_comp=obj.Fsquare_comp;
ldetC=obj.ldetC;
info=struct;

TStart=tic;
tStart=cputime;
optimality_found=0;

% create function call
func_alpha = @(x0,alpha) mix_DDFact_DDFact_comp_Linx_alpha(x0,C,s,...
    F,Fsquare,F_comp,Fsquare_comp,ldetC,...
    Gamma1,Gamma2,Gamma3,A_data,b_data,alpha(1),alpha(2),alpha(3));

%% check if any single upper bound is the best, i.e., if any unit alpha is optimal
% we leverage the facts that the mixing bound is convex in alpha

% check DDFact
% the DDFact bound and that the feasible direction at (1,0,0) can be
% written as the conic combination of (-1,1,0) and (-1,0,1), and so one so
% forth
alpha=[1,0,0]';
[fval,dalpha,info_mix] = func_alpha(x0,alpha);
% the critical point is that if the directional derivative is nonnegative
% in every feasible direction, then the current point is optimal. By the
% previous discussion, we only need to verify the nonnegativeness in two
% directions, (-1,1,0) and (-1,0,1)
if dalpha(2)-dalpha(1)>=0 && dalpha(3)-dalpha(1)>=0
    optimality_found=1;
    [~,~,info] = obj.DDFact_obj(info_mix.x,s,Gamma1);
    info.alpha=[1,0,0]';
end

% check DDFact comp
if optimality_found ==0
    alpha=[0,1,0]';
    [fval,dalpha,info_mix] = func_alpha(x0,alpha);
    if dalpha(1)-dalpha(2)>=0 && dalpha(3)-dalpha(2)>=0
        optimality_found=1;
        [~,~,info] = obj.DDFact_comp_obj(info_mix.x,s,Gamma2);
        info.alpha=[0,1,0]';
    end
end

% check Linx
if optimality_found ==0
    alpha=[0,0,1]';
    [fval,dalpha,info_mix] = func_alpha(x0,alpha);
    if dalpha(1)-dalpha(3)>=0 && dalpha(2)-dalpha(3)>=0
        optimality_found=1;
        [~,~,info] = obj.Linx_obj(info_mix.x,s,Gamma3);
        info.alpha=[0,0,1]';
    end
end

%% enter the search procedure that none of the single upper bound is optimal
if optimality_found ==0
    alpha=[1/3,1/3,1/3]'; % initialize alpha
    [fval,dalpha,info_mix] = func_alpha(x0,alpha);% evaluate objective value and gradient
    beta=0.8; % shrink parameter for t
    Tol=1e-10;
    Numiterations=100;
    k=0;
    % we use the proximal gradient method to optimize alpha over the unit
    % simplex @https://people.eecs.berkeley.edu/~elghaoui/Teaching/EE227A/lecture18.pdf
    nx=info_mix.x;
    while true
        k=k+1;
        %sprintf("iteration:%d, function value:%f",k,fval)
        t=1;
        while true % line search
            if t<1e-6
                break;
            end
            Gt=1/t*(alpha-proj_simplex_vector(alpha-t*dalpha));
            newalpha=alpha-t*Gt;
            [newfval,newdalpha,newinfo_mix] = func_alpha(nx,newalpha);
            if newfval<=fval-t*dalpha'*Gt+t/2*(Gt'*Gt)
                break
            end
            t=beta*t;
        end
        if k>Numiterations || (Gt'*Gt)<Tol ||fval-newfval<Tol
            alpha=newalpha;
            dalpha=newdalpha;
            fval=newfval;
            info_mix=newinfo_mix;
            break
        end
        alpha=newalpha;
        dalpha=newdalpha;
        fval=newfval;
        info_mix=newinfo_mix;
        nx=info_mix.x;
    end
    
    %% calculate dual solutions
    info.x=info_mix.x;
    info.dx=info_mix.dx;
    info.fval=fval;
    info.alpha=alpha;
    
    f=[zeros(n,1);ones(n,1);b_data;s];
    Aeq=[-eye(n),eye(n),A_data',ones(n,1)];
    beq=info.dx;
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

    info.continuous_dualgap=dualgap+info_mix.cache;
    info.dualbound=info.continuous_dualgap+info.fval;

    %% fixing variables
    info.fixnum=0;
    info.fixnum_to0=0;
    info.fixto0list=[];
    info.fixnum_to1=0;
    info.fixto1list=[];
    info.integrality_gap=info.dualbound-obj.obtain_lb(s);
    if info.integrality_gap>1e-6
        info.solved=0;
    else
        info.solved=1;
    end
    for i=1:n
        if info.integrality_gap<info.dual_upsilon(i)-1e-10
            info.fixnum=info.fixnum+1;
            info.fixnum_to0=info.fixnum_to0+1;
            info.fixto0list(end+1)=i;
        elseif info.integrality_gap<info.dual_nu(i)-1e-10
            info.fixnum=info.fixnum+1;
            info.fixnum_to1=info.fixnum_to1+1;
            info.fixto1list(end+1)=i;
        end
    end
end
fval=info.fval;
x=info_mix.x;
info.x=x;
time=toc(TStart);
tEnd=cputime-tStart;