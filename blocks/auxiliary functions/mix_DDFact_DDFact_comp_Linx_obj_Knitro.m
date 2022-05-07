function [fval,dx] = mix_DDFact_DDFact_comp_Linx_obj_Knitro(x,C,s,F,Fsquare,F_comp,Fsquare_comp,ldetC,Gamma1,Gamma2,Gamma3,alpha1,alpha2,alpha3)
[fval1,dx1] = DDFact_obj_Knitro(x,s,F,Fsquare,Gamma1);
[fval2,dx2] = DDFact_comp_obj_Knitro(x,s,F_comp,Fsquare_comp,ldetC,Gamma2);
[fval3,dx3] = Linx_obj_Knitro(x,C,Gamma3);
fval=alpha1*fval1+alpha2*fval2+alpha3*fval3;
dx=alpha1*dx1+alpha2*dx2+alpha3*dx3;

end

