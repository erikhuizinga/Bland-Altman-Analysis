function x = Ninv(p)
% Calculate inverse normal
x = -sqrt(2) .* erfcinv(2 * p);
end