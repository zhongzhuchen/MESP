%% obtain class properties
C=obj.C;
n = obj.size;
info=struct;
F = obj.F;
Fsquare = obj.Fsquare;
F_comp = obj.F_comp;
Fsquare_comp =obj.Fsquare_comp;

t1=tic;
if n>200
    TOL= 10^(-6);
    Numiterations=20; 
else
    TOL= 10^(-10);
    Numiterations=200; 
end

%% calculate the better lower bound among C and Cinv if C is invertible
% [U,D]=eig(C);
% lam=diag(D);
% if min(lam)> 0
%     shift=log(prod(lam));  % logdet C
%     Cinv=U*diag(1./lam)*U';
%     [~,heurval]=heur(C,n,s);        % HEURSITIC ON ORIGINAL
%     [~,cheurval]=heur(Cinv,n,n-s); % HEURISTIC ON COMPLEMENT
%     if cheurval+shift > heurval        % PICK THE BEST
%         heurval=cheurval+shift;
%     end
% else
%     [~,heurval]=heur(C,n,s);        % HEURSITIC ON ORIGINAL
% end
heurval = obj.obtain_lb(s);
x0=s/n*ones(n,1);
% Initialize Gamma
Gamma1=Gamma1Init;
Gamma2=Gamma2Init;

difgap=1;
k=1;

c1=1e-4;
c2=0.9;

%% calculate the gradient of mixing bound with respect to Gamma
if mix_pattern == "DDFact_Linx"
    [bound,x,info_mix]= obj.mix_DDFact_Linx(s,Gamma1,Gamma2);

    Fx=diag(sqrt(x))*F;
    Fsquarex=Fsquare.*reshape(x,1,1,n);
    
    [~,dGamma1,~] = DDFact_obj_auxiliary(Gamma1,s,Fx,Fsquarex);
    grad1=Gamma1.*dGamma1-x;

    scaleC=diag(Gamma2)*C;
    AUX = scaleC*diag(x)*scaleC';
    B = - diag(x) + eye(n);
    %Compute F(gamma,x)
    L= AUX+B;
    L=(L+L')/2; % force symmetry
    %Compute inv(F(gamma,x))
    [U,D]=eig(L);
    lam=diag(D);
    Finv=U*diag(1./lam)*U';
    %Compute the residual res
    grad2=diag(Finv*AUX)-x;
elseif mix_pattern == "DDFact_comp_Linx"
    [bound,x,info_mix]= obj.mix_DDFact_comp_Linx(s,Gamma1,Gamma2);

    y=ones(n,1)-x;
    Fy=diag(sqrt(y))*F;
    Fsquarey=Fsquare.*reshape(y,1,1,n);
    [~,dGamma1,~] = DDFact_obj_auxiliary(Gamma1,n-s,Fy,Fsquarey);
    grad1=Gamma1.*dGamma1-y;

    scaleC=diag(Gamma2)*C;
    AUX = scaleC*diag(x)*scaleC';
    B = - diag(x) + eye(n);
    %Compute F(gamma,x)
    L= AUX+B;
    L=(L+L')/2; % force symmetry
    %Compute inv(F(gamma,x))
    [U,D]=eig(L);
    lam=diag(D);
    Finv=U*diag(1./lam)*U';
    %Compute the residual res
    grad2=diag(Finv*AUX)-x;
elseif mix_pattern == "DDFact_DDFact_comp"
    [bound,x,info_mix]= obj.mix_DDFact_DDFact_comp(s,Gamma1,Gamma2);

    Fx=diag(sqrt(x))*F;
    Fsquarex=Fsquare.*reshape(x,1,1,n);
    [~,dGamma1,~] = DDFact_obj_auxiliary(Gamma1,s,Fx,Fsquarex);
    grad1=Gamma1.*dGamma1-x;

    y=ones(n,1)-x;
    Fy=diag(sqrt(y))*F;
    Fsquarey=Fsquare.*reshape(y,1,1,n);
    [~,dGamma2,~] = DDFact_obj_auxiliary(Gamma2,n-s,Fy,Fsquarey);
    grad2=Gamma2.*dGamma2-y;
else
    error("There is no such mixing pattern.")
end
alpha=info_mix.alpha;
grad=[(1-alpha)*grad1; alpha*grad2];
gap=bound-heurval;
res=norm(grad);

Gamma=[Gamma1;Gamma2];
allGamma=Gamma;
allres=res;
allbound=bound;

H=eye(2*n); % initialize the inverse Hessian approximation
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

    if mix_pattern == "DDFact_Linx"
        [nbound,nx,info_mix]= obj.mix_DDFact_Linx(s,Gamma1,Gamma2);

        Fx=diag(sqrt(nx))*F;
        Fsquarex=Fsquare.*reshape(nx,1,1,n);
        
        [~,dGamma1,~] = DDFact_obj_auxiliary(Gamma1,s,Fx,Fsquarex);
        grad1=Gamma1.*dGamma1-nx;
    
        scaleC=diag(Gamma2)*C;
        AUX = scaleC*diag(nx)*scaleC';
        B = - diag(nx) + eye(n);
        %Compute F(gamma,x)
        L= AUX+B;
        L=(L+L')/2; % force symmetry
        %Compute inv(F(gamma,x))
        [U,D]=eig(L);
        lam=diag(D);
        Finv=U*diag(1./lam)*U';
        %Compute the residual res
        grad2=diag(Finv*AUX)-nx;
    elseif mix_pattern == "DDFact_comp_Linx"
        [nbound,nx,info_mix]= obj.mix_DDFact_comp_Linx(s,Gamma1,Gamma2);

        y=ones(n,1)-nx;
        Fy=diag(sqrt(y))*F;
        Fsquarey=Fsquare.*reshape(y,1,1,n);
        [~,dGamma1,~] = DDFact_obj_auxiliary(Gamma1,n-s,Fy,Fsquarey);
        grad1=Gamma1.*dGamma1-y;
    
        scaleC=diag(Gamma2)*C;
        AUX = scaleC*diag(nx)*scaleC';
        B = - diag(nx) + eye(n);
        %Compute F(gamma,x)
        L= AUX+B;
        L=(L+L')/2; % force symmetry
        %Compute inv(F(gamma,x))
        [U,D]=eig(L);
        lam=diag(D);
        Finv=U*diag(1./lam)*U';
        %Compute the residual res
        grad2=diag(Finv*AUX)-nx;
    elseif mix_pattern == "DDFact_DDFact_comp"
        [nbound,nx,info_mix]= obj.mix_DDFact_DDFact_comp(s,Gamma1,Gamma2);
        Fx=diag(sqrt(x))*F;
        Fsquarex=Fsquare.*reshape(nx,1,1,n);
        [~,dGamma1,~] = DDFact_obj_auxiliary(Gamma1,s,Fx,Fsquarex);
        grad1=Gamma1.*dGamma1-nx;
    
        y=ones(n,1)-nx;
        Fy=diag(sqrt(y))*F;
        Fsquarey=Fsquare.*reshape(y,1,1,n);
        [~,dGamma2,~] = DDFact_obj_auxiliary(Gamma2,n-s,Fy,Fsquarey);
        grad2=Gamma2.*dGamma2-y;
    end
    alpha=info_mix.alpha;
    ngrad=[(1-alpha)*grad1; alpha*grad2];
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
        Gamma1=nGamma(1:n);
        Gamma2=nGamma((n+1):2*n);
        [nbound,x,info_mix]= mix_func(s,Gamma1,Gamma2);
        if mix_pattern == "DDFact_Linx"
            [nbound,nx,info_mix]= obj.mix_DDFact_Linx(s,Gamma1,Gamma2);
            Fx=diag(sqrt(nx))*F;
            Fsquarex=Fsquare.*reshape(nx,1,1,n);
            
            [~,dGamma1,~] = DDFact_obj_auxiliary(Gamma1,s,Fx,Fsquarex);
            grad1=Gamma1.*dGamma1-nx;
        
            scaleC=diag(Gamma2)*C;
            AUX = scaleC*diag(nx)*scaleC';
            B = - diag(nx) + eye(n);
            %Compute F(gamma,x)
            L= AUX+B;
            L=(L+L')/2; % force symmetry
            %Compute inv(F(gamma,x))
            [U,D]=eig(L);
            lam=diag(D);
            Finv=U*diag(1./lam)*U';
            %Compute the residual res
            grad2=diag(Finv*AUX)-nx;
        elseif mix_pattern == "DDFact_comp_Linx"
            [nbound,nx,info_mix]= obj.mix_DDFact_comp_Linx(s,Gamma1,Gamma2);
            y=ones(n,1)-nx;
            Fy=diag(sqrt(y))*F;
            Fsquarey=Fsquare.*reshape(y,1,1,n);
            [~,dGamma1,~] = DDFact_obj_auxiliary(Gamma1,n-s,Fy,Fsquarey);
            grad1=Gamma1.*dGamma1-y;
        
            scaleC=diag(Gamma2)*C;
            AUX = scaleC*diag(nx)*scaleC';
            B = - diag(nx) + eye(n);
            %Compute F(gamma,x)
            L= AUX+B;
            L=(L+L')/2; % force symmetry
            %Compute inv(F(gamma,x))
            [U,D]=eig(L);
            lam=diag(D);
            Finv=U*diag(1./lam)*U';
            %Compute the residual res
            grad2=diag(Finv*AUX)-nx;
        elseif mix_pattern == "DDFact_DDFact_comp"
            [nbound,nx,info_mix]= obj.mix_DDFact_DDFact_comp(s,Gamma1,Gamma2);
            
            Fx=diag(sqrt(x))*F;
            Fsquarex=Fsquare.*reshape(nx,1,1,n);
            [~,dGamma1,~] = DDFact_obj_auxiliary(Gamma1,s,Fx,Fsquarex);
            grad1=Gamma1.*dGamma1-nx;
        
            y=ones(n,1)-nx;
            Fy=diag(sqrt(y))*F;
            Fsquarey=Fsquare.*reshape(y,1,1,n);
            [~,dGamma2,~] = DDFact_obj_auxiliary(Gamma2,n-s,Fy,Fsquarey);
            grad2=Gamma2.*dGamma2-y;
        end
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
[bound1,~,~]= mix_func(s,ones(n,1),ones(n,1));
if optbound>bound1
    optbound=bound1;
    optGamma=ones(2*n,1);
end
info.optbound=optbound;
info.optGamma1=optGamma(1:n);
info.optGamma2=optGamma((n+1):(2*n));
Gamma1=info.optGamma1;
Gamma2=info.optGamma2;
time=toc(t1);
info.time=time;