load('data63.mat');
n=length(C);
s=10;
a=MESP(C,double.empty(0,n),double.empty(0,1));

% A = full(gallery('tridiag',n,1,1,1));
% b = 2*ones(n,1);