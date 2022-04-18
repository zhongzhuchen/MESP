%% obtain class properties
A_data=obj.A;
b_data=obj.b;
[m,~] = size(obj.A);
info=struct;
scaleC=diag(Gamma)*obj.C;
n = obj.size;

%% calling knitro
obj_fn =  @(x) obj.Linx_obj_knitro(x,s,Gamma);
lb=zeros(n,1);
ub=ones(n,1);
Aeq=ones(1,n);
beq=s;
A=A_data;
b=b_data;

if sum(abs(Aeq*x0-beq))>1e-10
    error('The initial point x0 is not feasible.')
end

TStart=tic;
tStart=cputime;

extendedFeatures.HessFcn = @(x,lambda) Linx_hessfun(x,lambda,scaleC);
options = knitro_options('algorithm', 3, 'convex', 1, 'derivcheck', 0, 'outlev', 0 , 'gradopt', 1, ...
                         'hessopt', 1, 'maxit', 1000, 'xtol', 1e-15, 'derivcheck_tol',1e-5,...
                         'feastol', 1e-10, 'opttol', 1e-10, 'bar_feasible',1,...
                         'bar_maxcrossit', 10);
[x,knitro_fval,exitflag,output,lambda,~] = knitro_nlp(obj_fn,x0,A,b,Aeq,beq,lb,ub,[],extendedFeatures,options);
% ========================================
time=toc(TStart);
tEnd=cputime-tStart;

%% assign values to info
info.exitflag=exitflag;
info.x=x; % optimal solution
info.knitro_fval = -knitro_fval;
[fval,dx,finalinfo] = obj.Linx_obj(x,s,Gamma);
info.dx=dx;
info.fval=fval;
info.continuous_dualgap=finalinfo.dualgap;
info.dualbound=finalinfo.dualbound;
info.dual_upsilon=finalinfo.dual_upsilon;
info.dual_nu=finalinfo.dual_nu;
info.dual_pi=finalinfo.dual_pi;
info.dual_tau=finalinfo.dual_tau;
info.ub_lambda=lambda.upper;
info.lb_lambda=lambda.lower;

info.time=time;
info.cputime=tEnd;
info.iterations=output.iterations;
info.funcCount=output.funcCount;
info.firstorderopt=output.firstorderopt;
info.constrviolation=output.constrviolation;
info.algorithm=output.algorithm;

%% fixing variables
info.integrality_gap=info.dualbound-obj.obtain_lb(s);
% number of fixing variables
info.fixnum=0;
info.fixnum_to0=0;
info.fixto0list=[];
info.fixnum_to1=0;
info.fixto1list=[];
intgap=info.integrality_gap;
if info.integrality_gap>1e-6
    info.solved=0;
else
    info.solved=1;
end
for i=1:n
    if intgap<info.dual_upsilon(i)-1e-10
        info.fixnum=info.fixnum+1;
        info.fixnum_to0=info.fixnum_to0+1;
        info.fixto0list(end+1)=i;
    elseif intgap<info.dual_nu(i)-1e-10
        info.fixnum=info.fixnum+1;
        info.fixnum_to1=info.fixnum_to1+1;
        info.fixto1list(end+1)=i;
    end
end

