% model1.M
function [fx,fxp,fy,fyp,fypyp,fypy,fypxp,fypx,fyyp,fyy,fyxp,fyx,fxpyp,fxpy,fxpxp,fxpx,fxyp,fxy,fxxp,fxx,f] = model1

syms DELTA ALFA BETTA G P LAMBDAZ ETA GAMA SIGM
syms c cp l lp k kp Z Zp

f1 = c + G*kp - (1-DELTA) * k - Z * k^ALFA * l^(1-ALFA);
f2 = c^(-SIGM) - (BETTA/G) * cp^(-SIGM) * ((1/P) * Zp * ALFA * kp^(ALFA-1)*lp^(1-ALFA) + 1 - DELTA);
f3 = GAMA * l^ETA * c^SIGM - (1/P) * (1-ALFA) * Z * k^ALFA * l^(-ALFA);
f4 = Zp - Z^LAMBDAZ;

f = [f1;f2;f3;f4];

x = [k Z];
y = [l c];
xp = [kp Zp];
yp = [lp cp];

% f = subs(f, [x,y,xp,yp], (exp([x,y,xp,yp])));

[fx,fxp,fy,fyp,fypyp,fypy,fypxp,fypx,fyyp,fyy,fyxp,fyx,fxpyp,fxpy,fxpxp,fxpx,fxyp,fxy,fxxp,fxx]=anal_deriv(f,x,y,xp,yp);