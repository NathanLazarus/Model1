function [KSTAR,CSTAR,LSTAR,WSTAR,RSTAR,GAMA,ETA]=model1_ss_numericsetGAMAandETA(DELTA,ALFA,BETTA,G,P,FRISCHELAS,STEADYSTATEL,SIGM,ZSTAR,u)
% This program computes the steady state 

sym_frisch_elas = FrischElasticity_symbolic(u);
sym_labor_supply = laborsupply(u);

cs = casadi.SX.sym('cs');
ls = casadi.SX.sym('ls');
ks = casadi.SX.sym('ks');
gama = casadi.SX.sym('gama');
eta = casadi.SX.sym('eta');

x = vertcat(cs, ls, ks, gama, eta);
x0 = ones([5 1])*0.5;
obj = 1;

nlp = struct('f', obj, 'x', x, 'g', constraint(cs,ls,ks,gama,eta,DELTA,ALFA,BETTA,G,P,FRISCHELAS,STEADYSTATEL,SIGM,ZSTAR,sym_frisch_elas,sym_labor_supply));

opts=struct;
opts.print_time=0;
opts.ipopt.print_level=0;
solver = casadi.nlpsol('solver', 'ipopt', nlp, opts);

sol = solver('x0', x0, 'lbg', -1e-8, 'ubg', 1e-8);

solution = full(sol.x(:,1));

CSTAR = solution(1);
LSTAR = solution(2);
KSTAR = solution(3);
GAMA = solution(4);
ETA = solution(5);
WSTAR=w_func(KSTAR,LSTAR,P,ZSTAR,ALFA);
RSTAR=little_r(KSTAR,LSTAR,P,ZSTAR,ALFA,DELTA);

end

function [constraintval] =  constraint(c,l,k,GAMA,ETA,DELTA,ALFA,BETTA,G,P,FRISCHELAS,STEADYSTATEL,SIGM,ZSTAR,sym_frisch_elas,sym_labor_supply)
 constraintval = ...
 	[c + G*k - (1-DELTA) * k - ZSTAR * k^ALFA * l^(1-ALFA);...
     1 - (BETTA/G) * big_R(k,l,P,ZSTAR,ALFA,DELTA);...
     eval(sym_labor_supply) + w_func(k,l,P,ZSTAR,ALFA);...
     l-STEADYSTATEL;...
     FRISCHELAS-eval(sym_frisch_elas)];
end

