#!/usr/bin/env wolframscript
(* ::Package:: *)

Remove["Global`*"];
cnow=c[k,\[Zeta],\[Epsilon]];
lnow = l[k,\[Zeta],\[Epsilon]];
kplus = (f[k,lnow,\[Zeta]]-cnow +(1- \[Delta])k)/G;
cplus = c[kplus,\[Zeta]plus,\[Epsilon]];
lplus = l[kplus,\[Zeta]plus,\[Epsilon]];
\[Zeta]plus=\[Zeta]^LambdaZ Exp[\[Epsilon] z];
bigRplus = \[Alpha] \[Zeta]plus kplus^(\[Alpha]-1) lplus^(1-\[Alpha])/P+1-\[Delta];
w = (1-\[Alpha])\[Zeta] k^\[Alpha] lnow^-\[Alpha]/P;
EulerEq = u'[cnow]-\[Beta] u'[G cplus]bigRplus;
LMkt=v'[lnow]-u'[cnow]w;
EQM={EulerEq,LMkt};
u[x_]=(1/(1-\[Sigma]))x^(1-\[Sigma]);v[y_]=(H/(1+\[Eta]))y^(1+\[Eta]);
f[x_,l_,\[Zeta]_] =l^(1-\[Alpha]) x^\[Alpha] \[Zeta];
ssk = 0.918117320596905; ssc = 0.342808192642141;H=21.449138149223;G=1.014;P=G^(1/(1-0.95));
ss={k->ssk,kplus->ssk,cplus -> ssc,lplus -> 0.3,\[Zeta]->1,\[Epsilon]->0,1.->1,0.->0,z^2->1,z->0};
functionargs = {k,\[Zeta],\[Epsilon]};
ssfunctionvals = functionargs//.ss;
c[Sequence@@ssfunctionvals]=ssc;
l[Sequence@@ssfunctionvals]=0.3;
\[Alpha]=0.32; LambdaZ=0.9;\[Sigma]=0.9;\[Beta]=0.98;\[Delta]=0.08;\[Eta]=2;
sol={};
derivs=D[EQM,k]//.ss;
vars = {Derivative[Sequence@@Join[{1},ConstantArray[0,Length[functionargs]-1]]][c][Sequence@@ssfunctionvals],Derivative[Sequence@@Join[{1},ConstantArray[0,Length[functionargs]-1]]][l][Sequence@@ssfunctionvals]};
PossibleSolutions = NSolve[derivs==0,vars];
D[kplus,k]//.ss/. 1.->1//. PossibleSolutions[[1]];
D[kplus,k]//.ss/. 1.->1//. PossibleSolutions[[2]];
solDk=Join[PossibleSolutions[[1]],{Derivative[Sequence@@Join[{1},ConstantArray[0,Length[functionargs]-1]]][kfunc][Sequence@@ssfunctionvals]-> D[kplus,k]//.ss//. PossibleSolutions[[1]]}];
order = 3;
unorderedTupleFunc[n_,m_]:=Subsets[Range[m+n-1],{n}]//Subtract[#,ConstantArray[Range[0,n-1],Length@#]]&
unorderedTuples = Flatten[Map[unorderedTupleFunc[#,Length[functionargs]]&,Range[order]],1];
asvariables = MapThread[#1-> #2&,{Range[Length[functionargs]],functionargs}];
outputs = unorderedTuples//.asvariables;
derivstotake = DeleteCases[outputs,{k}];
iszero = {k->0,\[Zeta]->0,\[Epsilon]->10};
checkforzeros = Map[Total,derivstotake//.iszero];
zeros = Pick[derivstotake,Inner[Or,Map[OddQ,Quotient[checkforzeros,10]],Map[OddQ,Quotient[checkforzeros,100]],list]];
myfunc[a_,b_]:=Join[a,Module[{derivs,vars},If[MemberQ[zeros,b],
listofpositions = Flatten[Map[Position[functionargs,#] &,b]];
whichderivs = ConstantArray[0,Length[functionargs]];
Map[whichderivs[[#]]+=1 &,listofpositions];
{Derivative[Sequence@@whichderivs][c][Sequence@@ssfunctionvals]->0,
Derivative[Sequence@@whichderivs][l][Sequence@@ssfunctionvals]->0,
Derivative[Sequence@@whichderivs][kfunc][Sequence@@ssfunctionvals]->0},
derivs=D[EQM,Sequence@@b]//.ss//.a;
listofpositions = Flatten[Map[Position[functionargs,#] &,b]];
whichderivs = ConstantArray[0,Length[functionargs]];
Map[whichderivs[[#]]+=1 &,listofpositions];
vars = {Derivative[Sequence@@whichderivs][c][Sequence@@ssfunctionvals],Derivative[Sequence@@whichderivs][l][Sequence@@ssfunctionvals]};
solution = NSolve[derivs==0,vars][[1]];
Join[solution,{Derivative[Sequence@@whichderivs][kfunc][Sequence@@ssfunctionvals]-> D[kplus,Sequence@@b]//.ss//.a//. solution}]]]]
NthOrderRules = Fold[myfunc,solDk,derivstotake];
vals = NthOrderRules[[All,2]];
primeIDs = MapThread[#1-> #2&,{functionargs,Map[Prime,Range[Length[functionargs]]]}];
outputIDs =Map[Times[Sequence@@ #]&,outputs//.primeIDs];
allpermutationsIDs = Map[Times[Sequence@@ #]&,Tuples[functionargs,order]//.primeIDs];
Flatten[Map[Position[outputIDs,#]&,allpermutationsIDs]];
towrite = Map[Flatten[Map[Position[outputIDs,#]&,Map[Times[Sequence@@ #]&,Tuples[functionargs,#]//.primeIDs]]]&,Range[order]];
klcmatrix = Map[{ Map[vals[[#]]&,3 #],Map[vals[[#]]&,3 #-1],Map[vals[[#]]&,3 #-2]}&,towrite];
Map[Export["C:/Users/Nathan/Downloads/PerturbationMethods/Model1/coefs"<>ToString[#]<>".csv",klcmatrix[[#]],"csv"]&,Range[order]];


