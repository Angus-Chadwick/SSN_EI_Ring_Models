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


N = 2;
delta = 0.0001;
lambda = -0.1;

for k=1:10000

c=0.02 * k;
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

N = 3;

uprime = ones([N,1]);
uprime = uprime / norm(uprime);
P = (1/sqrt(N)) * diag([1, (-lambda)^(-1), (-lambda)^(-2)]);
J = [lambda,1,0;0,lambda * 1.00,1;0,0,lambda * 1.00];
A = P * J * inv(P);


Sigma = lyap(A, Sigma_u);
rp = inv(A) * uprime;

FI = rp'* inv(Sigma) * rp;
FInorm = FI / ((-2/lambda)*uprime'* inv(Sigma_u) * uprime);
