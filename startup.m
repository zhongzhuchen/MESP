load('data2000.mat');
n=length(C);
m=5;
A1 = randi([1,10],m,n);
s=949;
x0=s/n*ones(n,1);
A = double.empty(0,n);
b = double.empty(0,1);
% [xind, ~] = heur(C,s,A,b);
% x=zeros(n,1);
% x(xind)=1;
% b1 = A1*x-ones(m,1);
Prob = MESP(C,A,b);
% gamma=Prob.BFGS_Linx_gamma(s);

% C=Prob.C;
% n = Prob.size;
% A_data=Prob.A;
% b_data=Prob.b;
% F=Prob.F;
% Fsquare=Prob.Fsquare;
% F_comp=Prob.F_comp;
% Fsquare_comp=Prob.Fsquare_comp;
% ldetC=Prob.ldetC;
% info=struct;
% 
% Gamma1=ones(n,1);
% Gamma2=ones(n,1);
% Gamma3=sqrt(gamma)*ones(n,1);
% 
% Prob1=MESP(C,A1,b1);
% 
% idx=[15
%     30
%     31
%     33
%     37
%     41
%     44
%     46
%     53
%     56
%     63];
% 
% A2=zeros(5,n);
% A2(1,idx(1:2))=1;
% A2(2,idx(3:4))=1;
% A2(3,idx(5:6))=1;
% A2(4,idx(7:8))=1;
% A2(5,idx(9:10))=1;
% 
% b2=ones(5,1);
% Prob2 = MESP(C,A2,b2);