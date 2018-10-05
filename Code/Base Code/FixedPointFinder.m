%% Find fixed point of network dynamics


%%

W = [network.connectivity.JEE, -network.connectivity.JEI; network.connectivity.JIE, -network.connectivity.JII];
u = [inputs.IE_FF', inputs.II_FF']';
T = diag([network.cells.tauE * ones(NE,1); network.cells.tauI * ones(NI,1)]);     
Tinv = inv(T);
gamma = 2;
r0 = R0';

% minimise function L 


Method = 'Newton';

if strcmp(Method, 'GradDesc')
    epsilon = 0.01;
      L = @(r, W, u, gamma) ( (eye(length(r)) - diag(gamma * max(0, W * r + u).^(gamma-1)) * W)' *  (r - max(0, W * r + u).^gamma) ); 
elseif strcmp(Method, 'Newton')
      L = @(r, W, u, gamma) ( (eye(length(r)) - diag(gamma * max(0, W * r + u).^(gamma-1)) * W) \  (r - max(0, W * r + u).^gamma) );  
epsilon = 1;
elseif strcmp(Method, 'NewtonSquared')
          L = @(s, W, u, gamma) ( 2 * (diag(s) - diag(gamma * s .* max(0, W * s.^2 + u).^(gamma-1)) * W) \  (s.^2 - max(0, W * s.^2 + u).^gamma) );  
epsilon = 1;
end



Niter = 100;
rmin = r0;

loss = 10;
i=1;

clear Ls

while and(loss > 1e-15, i<Niter)
    
    
    grad = L(rmin, W, u, gamma);
        
    rmin = rmin - epsilon * grad ;
 %   rmin = max(0,rmin);  % do we need this?
    loss = norm(Tinv * ( rmin - max(0, W * rmin + u).^2));  

 i = i+1;
 
 Ls(i-1) = loss;
    
end


 Niter_conv(nEI,nIE) = i;
 Niter_conv(nEI,nIE)

 
