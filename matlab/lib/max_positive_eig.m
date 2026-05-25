function lambda = max_positive_eig(matrix)
%MAX_POSITIVE_EIG Largest positive real eigenvalue.

vals = eig(matrix);
vals = real(vals(abs(imag(vals)) < 1e-7));
vals = vals(vals > 1e-9);
if isempty(vals)
    error('No positive real eigenvalue found.');
end
lambda = max(vals);
end
