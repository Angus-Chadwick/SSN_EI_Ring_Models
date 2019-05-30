%% Calculate 2d non-normal dynamical system Fisher Information

clear all

tau_a = 1;

qmax = 1000;

for q=1:qmax
    
    tau_b = tau_a * q/qmax * 10;

for i=1:21
    for j=1:21
        
theta_a = pi/2 * (i-11)/10;
theta_b = pi/2 * (j-11)/10;
theta_u = 0;

va = [cos(theta_a); sin(theta_a)];
vb = [cos(theta_b); sin(theta_b)];
uprime = [cos(theta_u); sin(theta_u)];

epsilon = va' * vb * 2 * sqrt(tau_a * tau_b) / (tau_a + tau_b);

I_F(i,j,q) = 1/(1-epsilon^2) * ((va'*uprime)^2 * tau_a * 2 + (vb'*uprime)^2 * tau_b * 2 - 4 * (epsilon * sqrt(tau_a * tau_b)) * (va'*uprime) * (vb'*uprime));

I_F_norm(i,j,q) = I_F(i,j,q) ./ ( 2 * max(tau_a,tau_b));

    end
end

end

for q=1:20
    subplot(5,4,q)
    
    s = q*10 - 9;

imagesc(90 * ([-10:10])/10, 90 * ([-10:10])/10, squeeze(I_F_norm(:,:,s)));
caxis([0,1.6])
colorbar
end



%%

clear all

tau_a = 1;

qmax = 1000;

Sigma_u = [1,0.5;0.5,1];

anglemax = 41;

for q=1:qmax
    
    tau_b = tau_a * q/qmax;

for i=1:anglemax
    for j=1:anglemax
       
theta_u = pi/2;
theta_a = pi/2 * (i-(anglemax+1)/2)/((anglemax-1)/2) + theta_u;
theta_b = pi/2 * (j-(anglemax+1)/2)/((anglemax-1)/2) + theta_u;

va = [cos(theta_a); sin(theta_a)];
vb = [cos(theta_b); sin(theta_b)];
uprime = [cos(theta_u); sin(theta_u)];


normchoice = 'Vleftnorm';

if strcmp(normchoice, 'Vleftnorm')

Vleft = [va'; vb'];
Vright = inv(Vleft);

elseif strcmp(normchoice, 'Vrightnorm')

Vright = [va, vb];
Vleft = inv(Vright);

end
Gamma = Vleft * Sigma_u * Vleft' ./ [2 * tau_a, tau_a + tau_b; tau_a + tau_b, 2 * tau_b]; 

FI(i,j,q) = uprime' * Vleft' * inv(Gamma) * Vleft * uprime;

end
end
end


FInorm = FI / (2*uprime'* inv(Sigma_u) * uprime);

for q=1:20
    subplot(5,4,q)
    
    s = (q-1) * (qmax/20) + 1;

imagesc(90 * ([-((anglemax-1)/2):((anglemax-1)/2)])/((anglemax-1)/2) + theta_u * 180/pi, 90 * ([-((anglemax-1)/2):((anglemax-1)/2)])/((anglemax-1)/2) + theta_u * 180/pi, squeeze(FInorm(:,:,s)));
title(strcat('\tau_b / \tau_a = ', num2str(s/qmax)))
colorbar
caxis([0,1.75])
xlabel('Angle (Mode 2)') ; ylabel('Angle (Mode 1)')
set(gca, 'ydir', 'normal')

end

subplot(5,4,11)
hold on
scatter(-5,5, 100, 'marker','*', 'linewidth', 2)

tau_b = 0.5;
theta_a = pi/2 * 5/180 + theta_u;
theta_b = pi/2 * (-5/180) + theta_u;

va = [cos(theta_a); sin(theta_a)];
vb = [cos(theta_b); sin(theta_b)];
uprime = [cos(theta_u); sin(theta_u)];

Vleft = [va'; vb'];

Gamma = Vleft * Sigma_u * Vleft' ./ [2 * tau_a, tau_a + tau_b; tau_a + tau_b, 2 * tau_b]; 
FI = uprime' * Vleft' * inv(Gamma) * Vleft * uprime;

J = inv(Vleft) * [-1/tau_a, 0; 0, -1/tau_b] * Vleft;
clear x
for i=1:100
t(i) = 0.1 * (i-1);
x1(:,i) = expm(J * 0.1 * (i-1)) * [0.1,0]';
x2(:,i) = expm(J * 0.1 * (i-1)) * [0.1,0.1]';
x3(:,i) = expm(J * 0.1 * (i-1)) * [0.0,0.1]';

end

figure 
subplot(2,3,1)
plot(t,x1')
xlabel('Time')
ylabel('Neuron Response')
legend('Neuron 1', 'Neuron 2')
subplot(2,3,2)
plot(t,x2')
xlabel('Time')
ylabel('Neuron Response')
legend('Neuron 1', 'Neuron 2')
subplot(2,3,3)
plot(t,x3')
xlabel('Time')
ylabel('Neuron Response')
legend('Neuron 1', 'Neuron 2')
subplot(2,3,4)
hold on
plot(t,va'*x1)
plot(t,vb'*x1)
xlabel('Time')
ylabel('Left Eigenvector Response')
legend('Evec 1', 'Evec 2')
subplot(2,3,5)
hold on
plot(t,va'*x2)
plot(t,vb'*x2)
xlabel('Time')
ylabel('Left Eigenvector Response')
legend('Evec 1', 'Evec 2')
subplot(2,3,6)
hold on
plot(t,va'*x3)
plot(t,vb'*x3)
xlabel('Time')
ylabel('Left Eigenvector Response')
legend('Evec 1', 'Evec 2')

%% 

for p=1:3

clearvars -except p

Sigma_u = [1,0;0,1];
theta_u = (p-1) * pi/4;

theta_offset = pi / 2 * 1/180 * 5;

tau_a = 1;
tau_b = 0.5;
theta_a = theta_offset + theta_u;
theta_b = -theta_offset + theta_u;

va = [cos(theta_a); sin(theta_a)];
vb = [cos(theta_b); sin(theta_b)];
uprime = [cos(theta_u); sin(theta_u)];
uprime_ortho = [-sin(theta_u); cos(theta_u)];

Vleft = [va'; vb'];

Gamma = Vleft * Sigma_u * Vleft' ./ [2 * tau_a, tau_a + tau_b; tau_a + tau_b, 2 * tau_b]; 
FI = uprime' * Vleft' * inv(Gamma) * Vleft * uprime;

J = inv(Vleft) * [-1/tau_a, 0; 0, -1/tau_b] * Vleft;


Vleftnorm = [uprime'; uprime_ortho'];
Jnorm = inv(Vleftnorm) * [-1/tau_a, 0; 0, -1/tau_b] * Vleftnorm;

noiseamp = 0.1;
tstep = 0.001;
Nt = 10000;
r1 = zeros([2,Nt]);
r2 = r1;
u1 = uprime_ortho * 0.01;
u2 = u1 + uprime * 0.01;
r1norm = zeros([2,Nt]);
r2norm = r1norm;

for i=1:(Nt-1)
    
    noise1 = [randn;randn];
    noise2 = [randn;randn];

    r1(:,i+1) = r1(:,i) + (J * r1(:,i) + u1 + noiseamp * noise1) * tstep;
    r2(:,i+1) = r2(:,i) + (J * r2(:,i) + u2 + noiseamp * noise2) * tstep;

    r1norm(:,i+1) = r1norm(:,i) + (Jnorm * r1norm(:,i) + u1 + noiseamp * noise1) * tstep;
    r2norm(:,i+1) = r2norm(:,i) + (Jnorm * r2norm(:,i) + u2 + noiseamp * noise2) * tstep;
    
end

Sigma_r = 0.5 * (cov(r1(:,(Nt/2):end)') + cov(r2(:,(Nt/2):end)'));
dmean = (mean(r2(:,(Nt/2):end)') - mean(r1(:,(Nt/2):end)'))';
rdisc = inv(Sigma_r) * dmean;
rdisc_ortho = [-dmean(2); dmean(1)];
rdisc = rdisc / norm(rdisc);
rdisc_ortho = rdisc_ortho / norm(rdisc_ortho);


subplot(4,3,p)
plot(r1')
hold on
plot(r2')
legend('Neuron 1 (stim 1)', 'Neuron 2 (stim 1)', 'Neuron 1 (stim 2)', 'Neuron 2 (stim 2)')
subplot(4,3,p + 3)
hold on

plot(rdisc' * r1)
plot(rdisc' * r2)
plot(rdisc_ortho' * r1)
plot(rdisc_ortho' * r2)

legend('Task-relevant Projection (stim 1)', 'Task-relevant Projection (stim 2)', 'Task-irrelevant Projection (stim 1)', 'Task-irrelevant projection (stim 2)')


Sigma_rnorm = 0.5 * (cov(r1norm(:,(Nt/2):end)') + cov(r2norm(:,(Nt/2):end)'));
dmeannorm = (mean(r2norm(:,(Nt/2):end)') - mean(r1norm(:,(Nt/2):end)'))';
rdiscnorm = inv(Sigma_rnorm) * dmeannorm;
rdisc_orthonorm = [-dmeannorm(2); dmeannorm(1)];
rdiscnorm = rdiscnorm / norm(rdiscnorm);
rdisc_orthonorm = rdisc_orthonorm / norm(rdisc_orthonorm);



subplot(4,3,p + 6)
plot(r1norm')
hold on
plot(r2norm')
legend('Neuron 1 (stim 1)', 'Neuron 2 (stim 1)', 'Neuron 1 (stim 2)', 'Neuron 2 (stim 2)')
subplot(4,3,p + 9)
hold on

plot(rdiscnorm' * r1norm)
plot(rdiscnorm' * r2norm)
plot(rdisc_orthonorm' * r1norm)
plot(rdisc_orthonorm' * r2norm)

legend('Task-relevant Projection (stim 1)', 'Task-relevant Projection (stim 2)', 'Task-irrelevant Projection (stim 1)', 'Task-irrelevant projection (stim 2)')


end