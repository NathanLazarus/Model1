function w = w_func(K,L,P,Z,ALFA)
w = (1./P) .* (1-ALFA) .* Z .* K.^ALFA .* L.^(-ALFA);
end