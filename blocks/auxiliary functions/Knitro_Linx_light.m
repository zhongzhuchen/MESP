function [fval,x,info] = Knitro_Linx_light(x0,C,s,A,b,Gamma)
%% obtain class properties
A_data=A;
b_data=b;
[m,~] = size(A);
scaleC=diag(Gamma)*C;
n = length(x0);
info = struct;
%% calling knitro
obj_fn =  @(x) Linx_obj_Knitro(x,C,Gamma);
lb=zeros(n,1);
ub=ones(n,1);
Aeq=ones(1,n);
beq=s;
A=A_data;
b=b_data;

if sum(abs(Aeq*x0-beq))>n*1e-10
    error('The initial point x0 is not feasible.')
end

TStart=tic;
tStart=cputime;

extendedFeatures.HessFcn = @(x,lambda) Linx_hessfun(x,lambda,scaleC);
options = knitro_options('algorithm', 4, 'convex', 1, 'derivcheck', 0, 'outlev', 0 , 'gradopt', 1, ...
                         'hessopt', 1, 'maxit', 1000, 'xtol', 1e-15, 'derivcheck_tol',1e-5,...
                         'feastol', 1e-10, 'opttol', 1e-10, 'bar_feasible',1,...
                         'bar_maxcrossit', 10);
[x,knitro_fval,exitflag,output,lambda,~] = knitro_nlp(obj_fn,x0,A,b,Aeq,beq,lb,ub,[],extendedFeatures,options);
% ========================================
time=toc(TStart);
tEnd=cputime-tStart;
fval=-knitro_fval;
info.fval=fval;
info.time=time;
info.cputime = tEnd;
end
