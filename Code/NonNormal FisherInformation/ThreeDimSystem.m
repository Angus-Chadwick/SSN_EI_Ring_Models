for k=1:1000
for j=1:1000
lambda = -1;
delta_b = 0.01;
delta_c = 0.02;
c_a = +0.02*(k-500);
c_b = +0.2*(j-500);
cba = +0.02*(k-500);
% c_a = +0.05*(k-500);
% c_b =  0.05*(j-500);
% cba = 0;


A = lambda * [1,-c_a, -c_b;0,1+delta_b,-cba;0,0,1+delta_c];

[V,D] = eig(A );
Vleft = inv(V);

Sigma_u = [1,0,0;0,1,0;0,0,1];

uprime = [0;1;1];
uprime = uprime / norm(uprime);

tau_a = -1/lambda;
tau_b = -1/(lambda * (1+delta_b));
tau_c = -1/(lambda * (1+delta_c));

Gamma = Vleft * Sigma_u * Vleft' ./ [2 * tau_a, tau_a + tau_b, tau_a + tau_c; tau_a + tau_b, 2 * tau_b, tau_b + tau_c; tau_c + tau_a, tau_c + tau_b, 2 * tau_c]; 
FI = uprime' * Vleft' * inv(Gamma) * Vleft * uprime;
FInorm(k,j) = FI / (2*uprime'* inv(Sigma_u) * uprime);

end
end

%%


clear FInorm


N = 3;
delta = 0.001;
lambda = -0.1;

for k=1:10000

c=0.01 * k;
%c = 1;


A = zeros(N);
A = A + diag([lambda:-delta:(lambda - (N-1) * delta)]) + c*(triu(ones(N),1) - triu(ones(N),2));
[V,D] = eig(A );
Vleft = inv(V);
Sigma_u = eye(N);
% uprime = zeros([N,1]);
% uprime(end) = 1;
 uprime = ones([N,1]);

uprime = uprime / norm(uprime);

tau = -1./diag(A);

Gamma = Vleft * Sigma_u * Vleft' ./ (tau + tau'); 

FI = uprime' * Vleft' * inv(Gamma) * Vleft * uprime;
FInorm(k) = FI / ((-2/lambda)*uprime'* inv(Sigma_u) * uprime);

end

for i=1:N
Vleft(i,:) = Vleft(i,:) ./ norm(Vleft(i,:));
end
Gamma = Vleft * Sigma_u * Vleft' ./ (tau + tau');
Delta = Gamma - mean(Gamma(:));
invDelta = inv(Delta);

%% Jordan normal form
for N=2:11

    clearvars -except N FInorm FInorm_wtest
    
lambda = -10;
Sigma_u = eye(N);

uprime = zeros([N,1]);
uprime(end) = 1;
uprime = uprime / norm(uprime);
c = 100000;
for i=1:N
    D(i) = c^(N-i);
end
P = diag(D);
J = eye(N) * lambda;
for i=1:(N-1)
    J(i,i+1) = 1;
end
A = P * J * inv(P);


Sigma = lyap(A, Sigma_u);
rp = inv(A) * uprime;

FI = rp'* inv(Sigma) * rp;
FInorm(N) = FI / ((-2/lambda)*uprime'* inv(Sigma_u) * uprime);

wout = inv(Sigma) * rp;

wtest = (- lambda / c).^fliplr(0:(N-1));

FI_wtest = (wtest * rp)^2 / (wtest * Sigma * wtest');
FInorm_wtest(N) = FI_wtest / ((-2/lambda)*uprime'* inv(Sigma_u) * uprime);  

end

% temporal integration using box filter approximation

for N=1:50

clear g1 g2

for j=1:N
for i=1:N
g1(i,j) = factorial(2*N - i - j) / (factorial(N - i) * factorial(N-j));
g2(i,j) = 2^-(2*N - i - j + 1);
end
end

g = g1.* g2;
TempInt(N) = N^2./sum(sum(g));

TempIntOpt(N) = sum(sum(inv(g)));

wopt{N}  = (inv(g) * ones([N,1])) .* (-lambda / c).^fliplr([0:(N-1)])';

end




%% cyclic network
clear all

tau = 10;

for N=3:16
    
    clear S

for s=1:1000

    clearvars -except N FInorm FIflatnorm FIoptnorm tau s c S 
    
lambda = -1/tau;
Sigma_u = eye(N);

uprime = zeros([N,1]);
uprime(end) = 1;
uprime = uprime / norm(uprime);
c(s) = 0.0005 * s;

D = c(s).^fliplr([0:(N-1)]);

% P = diag(D);
% J = eye(N) * lambda + diag(ones([N-1,1]), 1);
% A =P * J * inv(P);

A = lambda * eye(N) + diag(ones([N-1,1]), 1) * c(s);
A(N,1) = c(s);

try
Sigma = lyap(A, Sigma_u);
rp = inv(A) * uprime;

FI = rp'* inv(Sigma) * rp;
FInorm(s,N) = FI / ((-2/lambda)*uprime'* inv(Sigma_u) * uprime);
end


end
end
