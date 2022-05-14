function [fval,dalpha,info] = mix_DDFact_DDFact_comp_Linx_alpha(x0,C,s,...
    F,Fsquare,F_comp,Fsquare_comp,ldetC,...
    Gamma1,Gamma2,Gamma3,A_data,b_data,alpha1,alpha2,alpha3)
% Output the obj and gradient of the mixing bound at alpha
info=struct;
n=size(C,1);
obj_fn =  @(x) mix_DDFact_DDFact_comp_Linx_obj_Knitro(x,C,s,...
    F,Fsquare,F_comp,Fsquare_comp,ldetC,Gamma1,Gamma2,Gamma3,alpha1,alpha2,alpha3);
lb=zeros(n,1);
ub=ones(n,1);
Aeq=ones(1,n);
beq=s;
A=A_data;
b=b_data;

if sum(abs(Aeq*x0-beq))>n*1e-10
    error('The initial point x0 is not feasible.')
end

options = knitro_options('algorithm', 3, 'convex', 1, 'derivcheck', 0, 'outlev', 1 , 'gradopt', 1, ...
                         'hessopt', 2, 'maxit', 1000, 'xtol', 1e-15, ...
                         'feastol', 1e-10, 'opttol', 1e-10, 'bar_feasible',1,...
                         'bar_maxcrossit', 10);
TStart=tic;
tStart=cputime;
try
[x,knitro_fval,exitflag,output,lambda,~] = knitro_nlp(obj_fn,x0,A,b,Aeq,beq,lb,ub,[],[],options);  
catch
    1+1
end
time=toc(TStart);
tEnd=cputime-tStart;

%% objective value
fval=-knitro_fval;
%% gradient
[fval1, dx1, info1] = DDFact_obj_Knitro(x,s,F,Fsquare,Gamma1);
[fval2, dx2, info2] = DDFact_comp_obj_Knitro(x,s,F_comp,Fsquare_comp,ldetC,Gamma2);
[fval3, dx3, info3] = Linx_obj_Knitro(x,C,Gamma3);
dalpha=[-fval1;-fval2;-fval3];

info.fval=fval;
info.x=x;
info.dx=-alpha1*dx1-alpha2*dx2-alpha3*dx3;
info.cache = alpha1*(info1.cache+sum(x.*log(Gamma1)))...
    +alpha2*(info2.cache-sum(x.*log(Gamma2)))...
    +alpha3*(info3.cache+sum(x.*log(Gamma3)));
info.time=time;
info.cputime = tEnd;

end

