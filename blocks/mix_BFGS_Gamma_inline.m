%% obtain class properties
C=obj.C;
n = obj.size;
info=struct;
F = obj.F;
Fsquare = obj.Fsquare;
F_comp = obj.F_comp;
Fsquare_comp =obj.Fsquare_comp;
A_data = obj.A;
b_data = obj.b;
ldetC=obj.ldetC;

t1=tic;
if n>200
    TOL= 10^(-6);
    Numiterations=20; 
else
    TOL= 10^(-10);
    Numiterations=200; 
end

%% calculate the better lower bound among C and Cinv if C is invertible
gamma=gammaInit;
nx=s/n*ones(n,1);
heurval = obj.obtain_lb(s);

difgap=1;
k=1;

c1=1e-4;
c2=0.9;
%solve the linx ralaxation for gamma and obtain x
if mix_pattern == "DDFact_Linx"
    [bound,nx,info_mix]= mix_DDFact_Linx_light(nx,C,s,F,Fsquare,A_data,b_data,ones(n,1),sqrt(gamma)*ones(n,1));
elseif mix_pattern == "DDFact_comp_Linx"
    [bound,nx,info_mix]= mix_DDFact_comp_Linx_light(nx,C,s,F_comp,Fsquare_comp,ldetC,A_data,b_data,ones(n,1),sqrt(gamma)*ones(n,1));
end

AUX = C*diag(nx)*C;
%Compute F(gamma,x)
nF=gamma*AUX - diag(nx) + eye(n);
nF=(nF+nF')/2; % force symmetry
%Compute inv(F(gamma,x))
[U,D]=eig(nF);
lam=diag(D);
Finv=U*diag(1./lam)*U';
%Compute the residual res
gap=bound-heurval;
alpha=info_mix.alpha;
res= alpha*(n - s - trace(Finv*diag(ones(n,1)-nx))); 

allgamma=gamma;
allres=res;
allbound=bound;

while(k<=Numiterations && gap > TOL && abs(res) > TOL && difgap > TOL) 
    if k>1
        difgap=abs(allbound(k)-allbound(k-1));
    end
    %Compute MAT= H_G(gamma:=exp(psi))
    Mat = gamma*(ones(n,1)-nx)'*diag(Finv*AUX*Finv);
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
    if mix_pattern == "DDFact_Linx"
        [nbound,nx,info_mix]= mix_DDFact_Linx_light(nx,C,s,F,Fsquare,A_data,b_data,ones(n,1),sqrt(ngamma)*ones(n,1));
    elseif mix_pattern == "DDFact_comp_Linx"
        [nbound,nx,info_mix]= mix_DDFact_comp_Linx_light(nx,C,s,F_comp,Fsquare_comp,ldetC,A_data,b_data,ones(n,1),sqrt(ngamma)*ones(n,1));
    end
    nAUX = C*diag(nx)*C;
    nF=ngamma*nAUX - diag(nx) + eye(n);
    nF=(nF+nF')/2; % force symmetry
    [U,D]=eig(nF);
    lam=diag(D);
    nFinv=U*diag(1./lam)*U';
    nalpha=info_mix.alpha;
    nres= nalpha*(n - s - trace(nFinv*diag(ones(n,1)-nx))); 
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
        if mix_pattern == "DDFact_Linx"
            [nbound,nx,info_mix]= mix_DDFact_Linx_light(nx,C,s,F,Fsquare,A_data,b_data,ones(n,1),sqrt(ngamma)*ones(n,1));
        elseif mix_pattern == "DDFact_comp_Linx"
            [nbound,nx,info_mix]= mix_DDFact_comp_Linx_light(nx,C,s,F_comp,Fsquare_comp,ldetC,A_data,b_data,ones(n,1),sqrt(ngamma)*ones(n,1));
        end
        nAUX = C*diag(nx)*C;
        nF=ngamma*nAUX - diag(nx) + eye(n);
        nF=(nF+nF')/2; % force symmetry
        [U,D]=eig(nF);
        lam=diag(D);
        nFinv=U*diag(1./lam)*U';
        nalpha=info_mix.alpha;
        nres= nalpha*(n - s - trace(nFinv*diag(ones(n,1)-nx)));  
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
info.iterations=k-1;
info.gap=gap;
info.absres=abs(res);
info.difgap=difgap;
[optbound,optiteration]=min(allbound);
optgamma=allgamma(optiteration);
info.alpha=info_mix.alpha;
nx=s/n*ones(n,1);
if mix_pattern == "DDFact_Linx"
    [bound1,nx,info_mix]= mix_DDFact_Linx_light(nx,C,s,F,Fsquare,A_data,b_data,ones(n,1),ones(n,1));
elseif mix_pattern == "DDFact_comp_Linx"
    [bound1,nx,info_mix]= mix_DDFact_comp_Linx_light(nx,C,s,F_comp,Fsquare_comp,ldetC,A_data,b_data,ones(n,1),ones(n,1));
end
if optbound>bound1
    optbound=bound1;
    optgamma=1;
    info.alpha=info_mix.alpha;
end
info.optbound=optbound;
info.optgamma=optgamma;

time=toc(t1);
info.time=time;