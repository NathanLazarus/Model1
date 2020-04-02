function [KSTAR,CSTAR,LSTAR,WSTAR,RSTAR,GAMA]=model1_ss_numericsetGAMA(DELTA,ALFA,BETTA,G,P,ETA,STEADYSTATEL,SIGM,ZSTAR)
% This program computes the steady state 

cs = casadi.SX.sym('cstar', 1, 1);
ls = casadi.SX.sym('lstar', 1, 1);
ks = casadi.SX.sym('kstar', 1, 1);
gama = casadi.SX.sym('gama', 1, 1);


x = vertcat(cs, ls, ks,gama);
x0 = ones(4);
obj = 1;

nlp = struct('f', obj, 'x', x, 'g', constraint(cs,ls,ks,gama,DELTA,ALFA,BETTA,G,P,ETA,STEADYSTATEL,SIGM,ZSTAR));

opts=struct;
opts.print_time=0;
opts.ipopt.print_level=0;
S = casadi.nlpsol('S', 'ipopt', nlp, opts);

sol = S('x0', x0,...
            'lbg', -1e-8, 'ubg', 1e-8);

solution = full(sol.x(:,1));

CSTAR = solution(1);
LSTAR = solution(2);
KSTAR = solution(3);
GAMA = solution(4);
WSTAR=w_func(KSTAR,LSTAR,P,ZSTAR,ALFA);
RSTAR=little_r(KSTAR,LSTAR,P,ZSTAR,ALFA,DELTA);

end

function [constraintval] =  constraint(C,L,K,gama,DELTA,ALFA,BETTA,G,P,ETA,STEADYSTATEL,SIGM,ZSTAR)
 constraintval = [C + G*K - (1-DELTA) * K - ZSTAR * K^ALFA * L^(1-ALFA);...
     1 - (BETTA/G) * ((1/P) * ZSTAR * ALFA * K^(ALFA-1)*L^(1-ALFA) + 1 - DELTA);...
     gama * L^ETA * C^SIGM - (1/P) * (1-ALFA) * ZSTAR * K^ALFA * L^(-ALFA);...
     L-STEADYSTATEL];
end

