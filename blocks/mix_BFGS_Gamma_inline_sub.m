if mix_pattern == "DDFact_Linx"
    [nbound,nx,info_mix]= mix_DDFact_Linx_light(nx,s,Gamma1,Gamma2);

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
    [nbound,nx,info_mix]= mix_DDFact_comp_Linx_light(nx,s,Gamma1,Gamma2);

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
    [nbound,nx,info_mix]= mix_DDFact_DDFact_comp_light(nx,s,Gamma1,Gamma2);
    Fx=diag(sqrt(x))*F;
    Fsquarex=Fsquare.*reshape(nx,1,1,n);
    [~,dGamma1,~] = DDFact_obj_auxiliary(Gamma1,s,Fx,Fsquarex);
    grad1=Gamma1.*dGamma1-nx;

    y=ones(n,1)-nx;
    Fy=diag(sqrt(y))*F;
    Fsquarey=Fsquare.*reshape(y,1,1,n);
    [~,dGamma2,~] = DDFact_obj_auxiliary(Gamma2,n-s,Fy,Fsquarey);
    grad2=Gamma2.*dGamma2-y;
else
    error("There is no such a mixing pattern.")
end