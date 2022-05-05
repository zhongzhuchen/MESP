%% obtain class properties
C=obj.C;
n = obj.size;
info=struct;
F = obj.F;
Fsquare = obj.Fsquare;
F_comp = obj.F_comp;
Fsquare_comp =obj.Fsquare_comp;
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
heurval = obj.obtain_lb(s);
nx=s/n*ones(n,1);
% Initialize Gamma
Gamma1=Gamma1Init;
Gamma2=Gamma2Init;

difgap=1;
k=1;

c1=1e-4;
c2=0.9;

%% calculate the gradient of mixing bound with respect to Gamma
mix_BFGS_Gamma_inline_sub;
bound=nbound;
x=nx;

alpha=info_mix.alpha;

grad=[(1-alpha)*grad1; alpha*grad2];
gap=bound-heurval;
res=norm(grad);

Gamma=[Gamma1;Gamma2];
allGamma=Gamma;
allres=res;
allbound=bound;

% initialize the inverse Hessian approximation, force to be block diagonal
% in the following text
H=eye(2*n); 
%% we use the optimal solution of every last linx bound as the initial point 
% for solving the next linx bound, trick for accelarating optimization
nx=x;

%% loop
while(k<=Numiterations && gap > TOL && abs(res) > TOL && difgap > TOL)
    sprintf('iteration: %d',k)
    if k>1
        difgap=abs(allbound(k)-allbound(k-1));
        if k>=2 
            deltag1=deltag(1:n);
            deltag2=deltag((n+1):(2*n));
            deltay1=deltay(1:n);
            deltay2=deltay((n+1):(2*n));
            H1=H(1:n,1:n);
            H2=H((n+1):(2*n),(n+1):(2*n));
            s_gmuly1=deltag1'*deltay1;
            s_gmuly2=deltag2'*deltay2;
            if k==2
                % if alpha=0 or 1, we do not update the corresponding
                % Hessian approximation
                if alpha < 1-1e-6
                    gam1=s_gmuly1/(deltag1'*deltag1);
                     H1=gam1*H1;
                end
                if alpha >1e-6
                    gam2=s_gmuly2/(deltag2'*deltag2);
                    H2=gam2*H2;
                end
            end
            if(s_gmuly1>=TOL)
                Hg=H1*deltag1;
                Hgy=Hg*deltay1';
                gHg=deltag1'*Hg;
                yg=s_gmuly1;
                H1=H1+(yg+gHg)/(yg^2)*(deltay1*deltay1')-(Hgy+Hgy')/yg;
            end
            if(s_gmuly2>=TOL)
                Hg=H2*deltag2;
                Hgy=Hg*deltay2';
                gHg=deltag2'*Hg;
                yg=s_gmuly2;
                H2=H2+(yg+gHg)/(yg^2)*(deltay2*deltay2')-(Hgy+Hgy')/yg;
            end
            H=blkdiag(H1,H2);
        end
    end
    if rank(H)<2*n
        H=eye(2*n);
    end
    %Compute the search direction for Newton's method (dir=-res/Mat)
    dir=-H*grad;
    edir=exp(dir);
    if norm(edir)==Inf || norm(edir)==0 || isnan(norm(edir))
        if norm(dir)==Inf || isnan(norm(dir))
            dir=-grad;
        else
            dir=dir/norm(dir)*10;
        end   
    end
    %check if alfa=1 satisfies the Strong Wolfe Conditions
    alfa=1;
    nGamma=Gamma.*exp(alfa*dir);
    Gamma1=nGamma(1:n);
    Gamma2=nGamma((n+1):2*n);

    mix_BFGS_Gamma_inline_sub;

    alpha=info_mix.alpha;
    ngrad=[(1-alpha)*grad1; alpha*grad2];
    nres= norm(ngrad);

    %[U,D]=eig(H);
    % sprintf("nbound:%f, nres:%f, Hessmineig:%f, min(nGamma):%f", nbound, nres,min(diag(D)),min(nGamma))

    if nbound-bound>c1*alfa*dir'*grad
        judge=0;
    elseif abs(dir'*ngrad)>c2*abs(dir'*grad)
        judge=0;
    else
        judge=1;
    end
    %line search
    b=1;
    a=0;
    while judge==0
        alfa=(a+b)/2;
        nGamma=Gamma.*exp(alfa*dir);
        Gamma1=nGamma(1:n);
        Gamma2=nGamma((n+1):2*n);

        mix_BFGS_Gamma_inline_sub;

        alpha=info_mix.alpha;
        ngrad=[(1-alpha)*grad1; alpha*grad2];
        nres= norm(ngrad);
        if nbound-bound>c1*alfa*dir'*grad
            b=alfa;
        elseif abs(dir'*ngrad)>c2*abs(dir'*grad)
            a=alfa;
        else
            break
        end
        if abs(b-a)<1e-3
            break
        end   
    end
    deltay=log(nGamma./Gamma);
    Gamma=nGamma;
    deltag=ngrad-grad;
    grad=ngrad;
    res=nres;
    sprintf("res=%f", nres);
    x=nx;
    bound=nbound;
    gap=bound-heurval;
    allGamma=[allGamma,Gamma];
    allres=[allres,res];
    allbound=[allbound,bound];
    k=k+1; 
end
info.iterations=k-1;
info.alpha=info_mix.alpha;
info.gap=gap;
info.absres=abs(res);
info.difgap=difgap;
[optbound,optiteration]=min(allbound);
optGamma=allGamma(:,optiteration);
Gamma1=ones(n,1);
Gamma2=ones(n,1);
mix_BFGS_Gamma_inline_sub;
if optbound>nbound
    optbound=nbound;
    optGamma=ones(2*n,1);
end
info.optbound=optbound;
info.optGamma1=optGamma(1:n);
info.optGamma2=optGamma((n+1):(2*n));
Gamma1=info.optGamma1;
Gamma2=info.optGamma2;
time=toc(t1);
info.time=time;