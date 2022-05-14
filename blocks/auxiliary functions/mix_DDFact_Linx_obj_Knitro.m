function [fval,dx] = mix_DDFact_Linx_obj_Knitro(x,C,s,F,Fsquare,Gamma1,Gamma2,alpha)
% create a callback function for Knitro specifying objective value and gradient 
% This function calculate the objective value, gradient, and info
% of mixing Linx&DDFact bound
%{
Input:
x       - current point for the mix problem
C,s,F,Fsquare - the size of subset we want to choose, also equals to the summation of all elements of x
Gamma1  - symmetric diagonal scaling paramter for DDFact
Gamma2  - row diagonal scaling paramter for Linx
alpha   - mixing parameter

Output:
fval    - objective value at current point x
dx      - the gradient of the obejctive function at x
%}
try
if alpha>1-1e-8
    [fval,dx] =Linx_obj_Knitro(x,C,Gamma2);
elseif alpha<1e-8
    [fval,dx] = DDFact_obj_Knitro(x,s,F,Fsquare,Gamma1);
else
    [fval1,dx1] = DDFact_obj_Knitro(x,s,F,Fsquare,Gamma1);
    [fval2,dx2] = Linx_obj_Knitro(x,C,Gamma2);
    fval = (1-alpha)*fval1+ alpha*fval2;
    dx = (1-alpha)*dx1+ alpha*dx2;
end
catch
    pause;
end
end

