function [KSTAR,CSTAR,LSTAR,WSTAR,RSTAR]=model1_ss_numeric(DELTA,ALFA,BETTA,G,P,ETA,GAMA,SIGM,ZSTAR,u)
% This program computes the steady state 

sym_labor_supply = laborsupply(u);

cs = casadi.SX.sym('cs');
ls = casadi.SX.sym('ls');
ks = casadi.SX.sym('ks');


x = vertcat(cs, ls, ks);
x0 = ones([3 1]);
obj = 1;

nlp = struct('f', obj, 'x', x, 'g', constraint(cs,ls,ks,DELTA,ALFA,BETTA,G,P,ETA,GAMA,SIGM,ZSTAR,sym_labor_supply));

opts=struct;
opts.print_time=0;
opts.ipopt.print_level=0;
solver = casadi.nlpsol('solver', 'ipopt', nlp,opts);

sol = solver('x0', x0,'lbg', -1e-8, 'ubg', 1e-8);

solution = full(sol.x(:,1));

CSTAR = solution(1);
LSTAR = solution(2);
KSTAR = solution(3);
WSTAR = w_func(KSTAR,LSTAR,P,ZSTAR,ALFA);
RSTAR = little_r(KSTAR,LSTAR,P,ZSTAR,ALFA,DELTA);

end

function [constraintval] =  constraint(c,l,k,DELTA,ALFA,BETTA,G,P,ETA,GAMA,SIGM,ZSTAR,sym_labor_supply)
 constraintval = ...
 	[c + G*k - (1-DELTA) * k - ZSTAR * k^ALFA * l^(1-ALFA);...
     1 - (BETTA/G) * big_R(k,l,P,ZSTAR,ALFA,DELTA);...
     eval(sym_labor_supply) + w_func(k,l,P,ZSTAR,ALFA)];
end

