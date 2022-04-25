function [xind, heurval] = heur(C,s,A,b)
% a heuristic method for obtaining a good feasible solution to MESP
%{
Input:
C       - the data matrix C
s       - the size of the subset we want to choose
A,b     - parameters for the linear constraint: Ax <= b

Output:
xind    - indicies of the solution found
heurval - objective value at the solution found
%}
m = length(b);
n = length(C);
%% using gencontsol if there is not linear constraint
[contsol] = gencontsol_eig(C,s);
% greedily construct a discrete solution based on the continuous one
[~,sortind]=sort(contsol,"descend");
xind = sortind(1:s);
[xind,heurval]=localstep(C,xind); 
% greedily construct a discrete solution from scratch
[xindg,~] = greedy(C,s);
% greedily construct a discrete solution based on the greedy one
[indg,heurvalg]=localstep(C,xindg);
if heurvalg>heurval                    % choose the better bound
   heurval=heurvalg;
   xind=indg;
end

%% With linear constraint
if m>0
    xind1=xind;
    heurval1=heurval;

    [xind] = gencontsol_diag(C,s,A,b);
    [xind,heurval]=localstep_lin(C,xind,A,b);

    x = zeros(n,1);
    x(xind1)=1;
    if heurval1>heurval && sum(max(A*x-b,0))<=1e-12
       heurval=heurval1;
       xind=xind1;
    end
end

%% feasibility check
x = zeros(n,1);
x(xind)=1;
if abs(sum(x)-s)>1e-12
    error("Heuristic solution is infeasible.")
elseif m>0 && sum(max(A*x-b,0))>1e-12
    error("Heuristic solution is infeasible.")
end
end

