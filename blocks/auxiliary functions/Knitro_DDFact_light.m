function [fval,x,info]=Knitro_DDFact_light(x0,C,s,F,Fsquare,A,b,Gamma)
%% obtain class properties
n = length(x0);
A_data=A;
b_data=b;
info=struct;

%% calling knitro
obj_fn =  @(x) DDFact_obj_Knitro(x,s,F,Fsquare,Gamma);
lb=zeros(n,1);
ub=ones(n,1);
Aeq=ones(1,n);
beq=s;
A=A_data;
b=b_data;

if sum(abs(Aeq*x0-beq))>n*1e-10
    error('The initial point x0 is not feasible.')
end

options = knitro_options('algorithm', 0, 'convex', 1, 'derivcheck', 0, 'outlev', 0 , 'gradopt', 1, ...
                         'hessopt', 2, 'maxit', 1000, 'xtol', 1e-15, ...
                         'feastol', 1e-10, 'opttol', 1e-10, 'bar_feasible',1,...
                         'bar_maxcrossit', 10);
TStart=tic;
tStart=cputime;
[x,knitro_fval,exitflag,output,lambda,~] = knitro_nlp(obj_fn,x0,A,b,Aeq,beq,lb,ub,[],[],options);  
time=toc(TStart);
tEnd=cputime-tStart;

fval=-knitro_fval;
info.fval=fval;
info.time=time;
info.cputime = tEnd;
end
