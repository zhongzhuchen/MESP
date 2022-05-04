function [fval,dx] = mix_DDFact_DDFact_comp_obj_Knitro(x,s,F,Fsquare,F_comp,Fsquare_comp,ldetC,Gamma1,Gamma2,alpha)
% create a callback function for Knitro specifying objective value and gradient 
% This function calculate the objective value, gradient, and info
% of mixing DDFact&DDFact comp bound
%{
Input:
x       - current point for the mix problem
s       - the size of subset we want to choose, also equals to the summation of all elements of x
Gamma1  - symmtric diagonal scaling paramter for DDFact
Gamma2  - symmtric diagonal scaling paramter for DDFact comp
alpha   - mixing parameter

Output:
fval    - objective value at current point x
dx      - the gradient of the obejctive function at x
%}
[fval1,dx1] = DDFact_obj_Knitro(x,s,F,Fsquare,Gamma1);
[fval2,dx2] = DDFact_comp_obj_Knitro(x,s,F_comp,Fsquare_comp,ldetC,Gamma2);
fval = (1-alpha)*fval1+ alpha*fval2;
dx = (1-alpha)*dx1+ alpha*dx2;
end