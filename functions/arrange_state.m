function [phi, xphi, yphi, ps] = arrange_state(u, y)
%  find states of a nonlinear dynamic system to be approximated

phit = zeros(3, 1);

for t = 2 : length(y)
    tt = t - 1;
    xxt(:, tt) = [y(t-1), y(t), u(t)]';
end

[xx, ps] = mapminmax(xxt, 0, 1);
xx = xx';

phi = xx(1:end-1, :);
yphi = phi(:,1);

xphi = xx(2:end, 1:2);