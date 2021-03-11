function y = proddamas_support(D, x, S)

% fast product by Dtilde transpose D (product by the large DAMAS matrix)


z = D(:, S)*(x(S).* D(:, S)');
y = real(sum(conj(D) .* (z*D), 1).');


end