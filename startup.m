load('data124.mat');
n=length(C);
m=0;
A1 = randi([1,10],m,n);
s=55;
x0=s/n*ones(n,1);
A = double.empty(0,n);
b = double.empty(0,1);
[xind, ~] = heur(C,s,A,b);
x=zeros(n,1);
x(xind)=1;
b1 = A1*x-ones(m,1);
Prob = MESP(C,A,b);
gamma=Prob.BFGS_Linx_gamma(s);

C=Prob.C;
n = Prob.size;
A_data=Prob.A;
b_data=Prob.b;
F=Prob.F;
Fsquare=Prob.Fsquare;
F_comp=Prob.F_comp;
Fsquare_comp=Prob.Fsquare_comp;
ldetC=Prob.ldetC;
info=struct;

Gamma1=ones(n,1);
Gamma2=ones(n,1);
Gamma3=sqrt(gamma)*ones(n,1);