function var = decision_func_thirdorder(decision_matrix,statevars,ssvals,sigma_Z)

distancefromss = [statevars-ssvals,sigma_Z];
flatten = @(A) A(:);
var = decision_matrix*[1,distancefromss,flatten((distancefromss)' * (distancefromss))',flatten(distancefromss'*flatten((distancefromss)' * (distancefromss))')']';

end