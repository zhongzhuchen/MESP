function [xind,heurval]=heur(c,n,s)
  % CALCULATE A CONTINUOUS SOLUTION
  contsol=gencontsol(c,s);
  sum(contsol);
  % GREEDILY FIND A DISCRETE SOLUTION BASED ON THE CONTINUOUS ONE
  [sortcontsol,sortind]=sort(contsol);
  sortind=flipud(sortind);
  xind=sortind(1:s);       % initialize to continuous solution
  % xind=(1:s);            % or initialize to 1,...,s.
  % xind=randomind(n,s);   % or initialize randomly
  % initheurval=log(det(c(xind,xind)))  % initial value for heuristic
  %
  % RUN 2-OPT
  [xind,heurval]=localstep(c,xind);    
  %heurval                             % value after 2-opt

  % DO THE GREEDY ALGOITHM
  [indg,heurvalg]=greedy(c,n,s);
  %heurvalg
  % RUN 2-OPT
  [indg,heurvalg]=localstep(c,indg);    
  %heurvalg                           % value after 2-opt

% CHOOSE THE BETTER OF THESE TWO HEURISTICS
if heurvalg>heurval                    % choose the better bound
   heurval=heurvalg;
   xind=indg;
end
