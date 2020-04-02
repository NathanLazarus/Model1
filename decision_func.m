function var = decision_func(decision_matrix,statevars,ssvals,sigma_Z)

flatten = @(A) A(:);
var = decision_matrix*[1,statevars-ssvals,flatten((statevars-ssvals)' * (statevars-ssvals))',sigma_Z^2]';

end