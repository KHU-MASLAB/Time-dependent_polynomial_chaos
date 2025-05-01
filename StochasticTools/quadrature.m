function [loc,wts] = quadrature(dist,gauss_point)

if dist == "Uniform"
    gauss_opt = '5point';
    [wts, loc] = gaussQuadrature(gauss_opt); % gauss quadrature
elseif dist == "Normal"
    [loc,wts,~] = gengausshermquadrule2(gauss_point); % Normal dist.
    loc = loc'; wts = wts';
end