%% Find fixed point of network dynamics


%%

W = [network.connectivity.JEE, -network.connectivity.JEI; network.connectivity.JIE, -network.connectivity.JII];
u = [inputs.IE_FF', inputs.II_FF']';
T = diag([network.cells.tauE * ones(NE,1); network.cells.tauI * ones(NI,1)]);     
Tinv = inv(T);
gamma = network.cells.gamma;
r0 = R0';

% minimise function L 

if soft_thresh

    Method = 'SoftThresh';

else
    
    Method = 'NewtonEfficient';

end    
    
if strcmp(Method, 'GradDesc')
    epsilon = 0.01;
      L = @(r, W, u, gamma) ( (eye(length(r)) - diag(gamma * max(0, W * r + u).^(gamma-1)) * W)' *  (r - max(0, W * r + u).^gamma) ); 
elseif strcmp(Method, 'Newton')
      L = @(r, W, u, gamma) ( (eye(length(r)) - diag(gamma * max(0, W * r + u).^(gamma-1)) * W) \  (r - max(0, W * r + u).^gamma) );  
epsilon = 1;
elseif strcmp(Method, 'NewtonEfficient')
      L = @(r, W, X, gamma) ( (eye(length(r)) - diag(gamma * X.^(1-1/gamma)) * W) \  (r - X) );  % X = max(0, W * r + u).^gamma
      epsilon = 1;
elseif strcmp(Method, 'SoftThresh')
      L = @(r, W, X, gamma) ( (eye(length(r)) - diag(gamma * (1-exp(-X.^(1/gamma))) .* X.^(1-1/gamma)) * W) \  (r - X) );  % X = log(1+e^x)^gamma
      epsilon = 1;

end



Niter = 100; % normally 100
rmin = r0;

loss = 10;
losstol = 1e-15;
i=1;

clear Ls

if strcmp(Method, 'SoftThresh')

while and(loss > losstol, i<Niter)
    
     X = log(1+exp(W * rmin + u)).^gamma;

    grad = L(rmin, W, X, gamma);
        
    rmin = rmin - epsilon * grad ;
     loss = norm(Tinv * ( rmin - log(1+exp(W * rmin + u)).^gamma)); 

 i = i+1;
 
 Ls(i-1) = loss;
    
end

else 

    
while and(loss > losstol, i<Niter)
    
    X = max(0, W * rmin + u).^gamma;

    grad = L(rmin, W, X, gamma);
        
    rmin = rmin - epsilon * grad ;
    loss = norm(Tinv * ( rmin - max(0, W * rmin + u).^gamma));  

 i = i+1;
 
 Ls(i-1) = loss;
    
end

end

 Niter_conv(nEI,nIE) = i;
 
if StableSim(nEI,nIE) & Ls(i-1) > losstol
    
    disp('Failed to Find Fixed Point')
    
end
