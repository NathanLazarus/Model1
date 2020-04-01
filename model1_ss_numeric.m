function [KSTAR,CSTAR,LSTAR,WSTAR,RSTAR]=model1_ss_numeric(DELTA,ALFA,BETTA,G,P,ETA,GAMA,SIGM,ZSTAR)
% This program computes the steady state 


cs = casadi.SX.sym('cstar', 1, 1);
ls = casadi.SX.sym('lstar', 1, 1);
ks = casadi.SX.sym('kstar', 1, 1);


x = vertcat(cs, ls, ks);
x0 = ones(3);
obj = 1;

nlp = struct('f', obj, 'x', x, 'g', constraint(cs,ls,ks,DELTA,ALFA,BETTA,G,P,ETA,GAMA,SIGM,ZSTAR));

opts=struct;
opts.print_time=0;
opts.ipopt.print_level=0;
S = casadi.nlpsol('S', 'ipopt', nlp,opts);

sol = S('x0', x0,...
            'lbg', -1e-8, 'ubg', 1e-8);

solution = full(sol.x(:,1));

CSTAR = solution(1);
LSTAR = solution(2);
KSTAR = solution(3);
WSTAR=ZSTAR/P*(1-ALFA)*KSTAR^ALFA*LSTAR^(-ALFA);
RSTAR=ZSTAR/P*ALFA*KSTAR^(ALFA-1)*LSTAR^(1-ALFA)-DELTA;

end

function [constraintval] =  constraint(C,L,K,DELTA,ALFA,BETTA,G,P,ETA,GAMA,SIGM,ZSTAR)
 constraintval = [C + G*K - (1-DELTA) * K - ZSTAR * K^ALFA * L^(1-ALFA);...
     1 - (BETTA/G) * ((1/P) * ZSTAR * ALFA * K^(ALFA-1)*L^(1-ALFA) + 1 - DELTA);...
     GAMA * L^ETA * C^SIGM - (1/P) * (1-ALFA) * ZSTAR * K^ALFA * L^(-ALFA)];
end

