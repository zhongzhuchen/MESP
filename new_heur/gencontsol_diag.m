function [xind] = gencontsol_diag(C,s,A,b)
n=size(C,1);
f=log(diag(C));
Aeq=ones(1,n);
beq=s;
xType=2*ones(n,1);
lb=zeros(n,1);
ub=ones(n,1);
x0=[];
options = knitro_options('algorithm',3,...  % active-set/simplex algorithm
                         'outlev',0);       % iteration display
[x,heurval,exitflag,~] = ...
    knitro_milp(f,xType,A,b,Aeq,beq,lb,ub,x0,[],options);
if exitflag<-199
    error('The MESP instance is infeasible.')
end
xind=transpose(1:n);
xind=xind(x>0);
end

