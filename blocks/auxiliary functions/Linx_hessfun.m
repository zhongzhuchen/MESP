function [H]= Linx_hessfun(x,lambda,scaleC)
    n=length(x);
    F=scaleC*diag(x)*scaleC'+eye(n)-diag(x);
    F=F+F';
    Finv=inv(F);
    hess1=-Finv.^2;
    hess21=Finv*scaleC; 
    hess22=hess21.^2;
    hess2=hess22+hess22';
    hess31=scaleC'*hess21;
    hess3=-hess31.^2;
    H=-(hess3+hess2+hess1);
    H=(H+H');
end

