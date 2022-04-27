%% obtain class properties
C=obj.C;
n = obj.size;
info=struct;

t1=tic;
%% setting different numbers of iterations for different problem size
if n>200
    TOL= 10^(-6);
    Numiterations=20; 
else
    TOL= 10^(-10);
    Numiterations=200; 
end

%% calculate the better lower bound among C and Cinv if C is invertible
heurval = obj.obtain_lb(s);

F = obj.F_comp;
Fsquare = obj.Fsquare_comp;
x0=s/n*ones(n,1);
%% Initializes Gamma and hyper parameter
Gamma=GammaInit;

difgap=1;
k=1;

c1=1e-4;
c2=0.9;

%% calculate the gradient of fact bound with respect to Gamma
[bound,x,~]= obj.Knitro_DDFact_comp(x0,s,Gamma);

y=ones(n,1)-x;
Fy=diag(sqrt(y))*F;
Fsquarey=Fsquare.*reshape(y,1,1,n);

[~,dGamma,~] = DDFact_obj_auxiliary(Gamma,n-s,Fy,Fsquarey);

gap=bound-heurval;
grad=Gamma.*dGamma-y;
res= norm(grad);

allGamma=Gamma;
allres=res;
allbound=bound;

H=eye(n); % initialize the inverse Hessian approximation
% sprintf('k: %d, gap: %f, abs(res): %f, difgap: %f',k, gap, abs(res), difgap)

%% we use the optimal solution of every last linx bound as the initial point 
% for solving the next linx bound, trick for accelarating optimization
nx=x;

%% loop
while(k<=Numiterations && gap > TOL && abs(res) > TOL && difgap > TOL)
    % sprintf('iteration: %d',k);
    if k>1
        difgap=abs(allbound(k)-allbound(k-1));
        if k>=2 
            s_gmuly=deltag'*deltay;
            if k==2 
                gam=s_gmuly/(deltag'*deltag);
                H=gam*H;
            end
            if(s_gmuly>=TOL)
                Hg=H*deltag;
                Hgy=Hg*deltay';
                gHg=deltag'*Hg;
                yg=s_gmuly;
                H=H+(yg+gHg)/(yg^2)*(deltay*deltay')-(Hgy+Hgy')/yg;
            end
        end
    end
    if rank(H)<n
        H=eye(n);
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
    [nbound,nx,~]= obj.Knitro_DDFact_comp(nx,s,nGamma);
    ny=ones(n,1)-nx;
    nFy=diag(sqrt(ny))*F;
    nFsquarey=Fsquare.*reshape(ny,1,1,n);
    
    [~,dnGamma,~] = DDFact_obj_auxiliary(nGamma,n-s,nFy,nFsquarey);

    ngrad=nGamma.*dnGamma-ny;
    nres= norm(ngrad);
    [U,D]=eig(H);
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
        [nbound,nx,~]= obj.Knitro_DDFact_comp(nx,s,nGamma);
        ny=ones(n,1)-nx;
        nFy=diag(sqrt(ny))*F;
        nFsquarey=Fsquare.*reshape(ny,1,1,n);
        
        [~,dnGamma,~] = DDFact_obj_auxiliary(nGamma,n-s,nFy,nFsquarey);

        ngrad=nGamma.*dnGamma-ny;
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
    x=nx;
    bound=nbound;
    gap=bound-heurval;
    allGamma=[allGamma,Gamma];
    allres=[allres,res];
    allbound=[allbound,bound];
    k=k+1; 

    %sprintf('k: %d, gap: %f, abs(res): %f, difgap: %f',k, gap, abs(res), difgap)
end
info.iterations=k-1;
info.gap=gap;
info.absres=abs(res);
info.difgap=difgap;
[optbound,optiteration]=min(allbound);
optGamma=allGamma(:,optiteration);
[bound1,~,~]= obj.Knitro_DDFact_comp(x0,s,ones(n,1));
if optbound>bound1
    optbound=bound1;
    optGamma=ones(n,1);
end
info.optbound=optbound;
info.optGamma=optGamma;
time=toc(t1);
info.time=time;