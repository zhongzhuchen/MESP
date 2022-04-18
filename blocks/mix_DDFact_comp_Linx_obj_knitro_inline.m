[fval1,dx1,~] = obj.DDFact_comp_obj(x,s,Gamma1);
[fval2,dx2,~] = obj.Linx_obj(x,s,Gamma2);
fval = -((1-alpha)*fval1+ alpha*fval2);
dx = -((1-alpha)*dx1+ alpha*dx2);