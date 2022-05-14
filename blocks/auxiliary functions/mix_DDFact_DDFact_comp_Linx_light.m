function  [fval,x,info] = mix_DDFact_DDFact_comp_Linx_light(x0,C,s,...
    F,Fsquare,F_comp,Fsquare_comp,ldetC,...
    Gamma1,Gamma2,Gamma3,A_data,b_data)
%% obtain class properties
n = length(x0);
m = size(A_data,1);
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
    info.alpha=[1,0,0]';
end

% check DDFact comp
if optimality_found ==0
    alpha=[0,1,0]';
    [fval,dalpha,info_mix] = func_alpha(x0,alpha);
    if dalpha(1)-dalpha(2)>=0 && dalpha(3)-dalpha(2)>=0
        optimality_found=1;
        info.alpha=[0,1,0]';
    end
end

% check Linx
if optimality_found ==0
    alpha=[0,0,1]';
    [fval,dalpha,info_mix] = func_alpha(x0,alpha);
    if dalpha(1)-dalpha(3)>=0 && dalpha(2)-dalpha(3)>=0
        optimality_found=1;
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
            [newfval,newdalpha,newinfo_mix] = func_alpha(x0,newalpha);
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
    
end
fval=info_mix.fval;
x=info_mix.x;
info.x=x;
time=toc(TStart);
tEnd=cputime-tStart;
info.time=time;
info.cputime=tEnd;
end

