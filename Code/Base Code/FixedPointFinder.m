%% Find fixed point of network dynamics

% inputs 
% r0 - initial guess at FP
% W - weights
% u - inputs
% gamma - supralinear exponent

W = [network.connectivity.JEE, -network.connectivity.JEI; network.connectivity.JIE, -network.connectivity.JII];
u = [inputs.IE_FF', inputs.II_FF']';
gamma = 2;
r0 = R0';

% minimise function L (try norm of L first)

L = @(r, W, u, gamma) ( (r - max(0, W * r + u).^gamma)' *  (eye(length(r)) - diag(gamma * max(0, W * r + u).^(gamma-1)) * W) );  % need to add a lagrange multiplier

epsilon = 0.01;
rmin = r0;

for i=1:1000
    
    rmin = rmin - epsilon * L(rmin, W, u, gamma)';
    rmin = max(0,rmin);
    %loss(i) = norm( L(rmin, W, u, gamma)); 
    loss(i) = norm(inv(T) * ( rmin - max(0, W * rmin + u).^2));
    loss(i)
    
    %epsilon = epsilon * 0.98;

end
