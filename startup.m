load('data63.mat');
n=length(C);
m=5;
% A1 = randi([1,10],m,n);
% s=940;
A = double.empty(0,n);
b = double.empty(0,1);
% [xind, ~] = heur(C,s,A,b);
% x=zeros(n,1);
% x(xind)=1;
% b1 = A1*x-ones(m,1);
Prob = MESP(C,A,b);
