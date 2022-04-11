function [F,Fsquare,ldetC] = gen_data(C,comp)
% generate data needed for calculating the objective value of DDFact
%{
Input:
C         - data matrix
comp      - equals to 1 means we caculate the complementary bound,
            otherwise caculate the original bound

Output:
C         - data matrix n-by-n
F         - a factor of C such that C=F*F', n-by-d
Fsquare   - a 3d array where Fsquare(:,:,i) represents the F(i,:)'*F(i,:)
ldetC     - if comp=1, ldetC=log det(inv(C)); if comp=0, ldetC=log det(C)
            and will not be used
%}
n=length(C);
ldetC=log(det(C));
if comp==1
    C=inv(C);
end
d=rank(C);
[U,D]=eig(C);
s=diag(D);
if ~issorted(s,'descend')
    [s,I] = sort(s, 'descend');
    U = U(:, I);
F=diag(sqrt(s))*U';
F=F(1:d,1:n);
F=F';
Fsquare=zeros(d,d,n);

for i=1:n
    Fsquare(:,:,i)=F(i,:)'*F(i,:);
end
end

