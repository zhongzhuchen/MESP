%% obtain class properties
C=obj.C;
n = obj.size;
A_data=obj.A;
b_data=obj.b;
[m,~] = size(obj.A);
info=struct;

%% solve linx and fact bound
x0=s/n*ones(n,1);
alpha=0;
% search interval
a=0;
b=1;
TStart=tic;
tStart=cputime;
% obtain optimal solution at search boundary a=0, b=1
[~, xa, ~] = obj.Knitro_DDFact(x0,s,Gamma1);
[~, xb, ~] = obj.Knitro_Linx(x0,s,Gamma2);
[fval_a1, dx_a1, info_a1] = obj.DDFact_obj(xa,s,Gamma1);
[fval_b2, dx_b2, info_b2] = obj.Linx_obj(xb,s,Gamma2);
[fval_a2, dx_a2, info_a2] = obj.Linx_obj(xa,s,Gamma2);
[fval_b1, dx_b1, info_b1] = obj.DDFact_obj(xb,s,Gamma1);
% gradient of the mixing objective function
dx_a = dx_a1;
dx_b = dx_b2;
% gradient of the mixing parameter alpha
dalpha_a = fval_a2-fval_a1;
dalpha_b = fval_b2-fval_b1;

if dalpha_a>=0
    fval = fval_a1;
    x = xa;
    dx = dx_a;
    alpha = 0;
    info1 = info_a1;
    info2 = info_a2;
elseif dalpha_b <=0
    fval = fval_b2;
    x = xb;
    dx = dx_b;
    alpha = 1;
    info1 = info_b1;
    info2 = info_b2;
else
    lb=zeros(n,1);
    ub=ones(n,1);
    Aeq=ones(1,n);
    beq=s;
    options = knitro_options('algorithm', 3, 'convex', 1, 'derivcheck', 0, 'outlev', 0 , 'gradopt', 1, ...
                             'hessopt', 2, 'maxit', 1000, 'xtol', 1e-15, ...
                             'feastol', 1e-10, 'opttol', 1e-10, 'bar_feasible',1,...
                             'bar_maxcrossit', 10);
    while b-a > 1e-6
        x0=(xa+xb)/2;
        if sum(abs(Aeq*x0-beq))>1e-10
            error('The initial point x0 is not feasible.')
        end
        alpha=(a+b)/2;
        obj_fn = @(x) obj.mix_DDFact_Linx_obj_knitro(x,s,Gamma1,Gamma2,alpha);
        [x,knitro_fval,exitflag,output,lambda,~] = knitro_nlp(obj_fn,x0,A_data,b_data,Aeq,beq,lb,ub,[],[],options);
        knitro_fval = - knitro_fval;
        [fval1, dx1, info1] = obj.DDFact_obj(x,s,Gamma1);
        [fval2, dx2, info2] = obj.Linx_obj(x,s,Gamma2);
        fval = (1-alpha)*fval1+ alpha*fval2;
        dx = (1-alpha)*dx1+(alpha)*dx2;
        dalpha = fval2 - fval1;
        if abs(dalpha)<1e-8
            break
        elseif dalpha<0
            a=alpha;
            xa=x;
        else
            b=alpha;
            xb=x;
        end
    end
end
time=toc(TStart);
tEnd=cputime-tStart;

%% assign values to info
% record important information
info.x=x; % optimal solution
info.dx=dx;
info.fval=fval;
info.alpha = alpha;

% calculate dual solution
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

% calculate continuous dualgap
info.continuous_dualgap=dualgap+(1-alpha)*(info1.cache+sum(x.*log(Gamma1)))+alpha*(info2.cache+sum(x.*log(Gamma2)));
info.dualbound=info.continuous_dualgap+fval;
info.time=time;
info.cputime=tEnd;

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
