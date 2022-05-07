%% obtain class properties
C=obj.C;
n = obj.size;
A_data=obj.A;
b_data=obj.b;
[m,~] = size(obj.A);
F=obj.F;
Fsquare=obj.Fsquare;
F_comp=obj.F_comp;
Fsquare_comp=obj.Fsquare_comp;
ldetC=obj.ldetC;
info=struct;

TStart=tic;
tStart=cputime;

if n>200
    TOL= 10^(-6);
    Numiterations=20; 
else
    TOL= 10^(-10);
    Numiterations=100; 
end


% Naive version 
func_Gamma=@(x0,gamma) mix_DDFact_DDFact_comp_Linx_light(x0,C,s,...
    F,Fsquare,F_comp,Fsquare_comp,ldetC,...
    ones(n,1),ones(n,1),sqrt(gamma)*ones(n,1),A_data,b_data);
%% calculate the better lower bound among C and Cinv if C is invertible
heurval = obj.obtain_lb(s);
nx=s/n*ones(n,1);
x0=nx;
% Initialize Gamma
gamma=gammaInit;
ngamma=gamma;
difgap=1;
k=1;

c1=1e-4;
c2=0.9;

[nbound,nx,info_mix]= func_Gamma(nx,ngamma);
alpha=info_mix.alpha;

AUX = C*diag(nx)*C;
%Compute F(gamma,x)
F=gamma*AUX - diag(nx) + eye(n);
F=(F+F')/2; % force symmetry
%Compute inv(F(gamma,x))
[U,D]=eig(F);
lam=diag(D);
Finv=U*diag(1./lam)*U';
%Compute the residual res
bound=nbound;
x=nx;

gap=bound-heurval;
res= alpha(3)*(n - s - trace(Finv*diag(ones(n,1)-x))); 
allgamma=gamma;
allres=res;
allbound=bound;

nx=x;
while(k<=Numiterations && gap > TOL && abs(res) > TOL && difgap > TOL) 
    if k>1
        difgap=abs(allbound(k)-allbound(k-1));
    end
    %Compute MAT= H_G(gamma:=exp(psi))
    Mat = gamma*(ones(n,1)-x)'*diag(Finv*AUX*Finv);
    %Compute the search direction for Newton's method (dir=-res/Mat)
    dir=-res/Mat;
    edir=exp(dir);
    if abs(edir)==Inf
        edir=sign(edir)*1e10;
    end
    dir=log(edir);
    %check if alfa=1 satisfies the Strong Wolfe Conditions
    alfa=1;
    ngamma=gamma*exp(alfa*dir);
    [nbound,nx,info_mix]= func_Gamma(nx,ngamma);
    alpha=info_mix.alpha;
    nAUX = C*diag(nx)*C;
    nF=ngamma*nAUX - diag(nx) + eye(n);
    nF=(nF+nF')/2; % force symmetry
    [U,D]=eig(nF);
    lam=diag(D);
    nFinv=U*diag(1./lam)*U';
    nres= alpha(3)*(n - s - trace(nFinv*diag(ones(n,1)-nx))); 
    if nbound-bound>c1*alfa*dir*res
        judge=0;
    elseif abs(dir*nres)>c2*abs(dir*res)
        judge=0;
    else
        judge=1;
    end
    %line search
    b=1;
    a=0;
    while judge==0
        alfa=(a+b)/2;
        ngamma=gamma*exp(alfa*dir);
        [nbound,nx,info_mix]= func_Gamma(nx,gamma);
        alpha=info_mix.alpha;
        nAUX = C*diag(nx)*C;
        nF=ngamma*nAUX - diag(nx) + eye(n);
        nF=(nF+nF')/2; % force symmetry
        [U,D]=eig(nF);
        lam=diag(D);
        nFinv=U*diag(1./lam)*U';
        nres= alpha(3)*(n - s - trace(nFinv*diag(ones(n,1)-nx))); 
        if nbound-bound>c1*alfa*dir*res
            b=alfa;
        elseif abs(dir*nres)>c2*abs(dir*res)
            a=alfa;
        else
            break
        end
        if abs(b-a)<1e-3
            break
        end   
    end
    gamma=ngamma;
    res=nres;
    Finv=nFinv;
    AUX=nAUX;
    x=nx;
    bound=nbound;
    gap=bound-heurval;
    allgamma=[allgamma,gamma];
    allres=[allres,res];
    allbound=[allbound,bound];
    k=k+1; 
end
[optbound,optiteration]=min(allbound);
optgamma=allgamma(optiteration);
[bound1,~,~] = func_Gamma(x0,1);
if optbound>bound1
    optbound=bound1;
    optgamma=1;
end
[fval,x,info]=obj.mix_DDFact_DDFact_comp_Linx(x0,s,ones(n,1),ones(n,1),sqrt(optgamma)*ones(n,1));
info.optgamma=optgamma;
info.iterations=k-1;
info.absres=abs(res);
info.difgap=difgap;
time=toc(TStart);
tEnd=cputime-tStart;
info.time=time;
info.cputime=tEnd;