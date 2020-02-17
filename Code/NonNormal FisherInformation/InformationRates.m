%% get rate of information gain for linear projection over different filters

clear all

lambda = - 1;
tau = -1/lambda;

t = 0:0.01:25;

%% signal processing filters

% box filter

Fbox = heaviside(tau-t); % filter

Ibox = (t - tau) .* heaviside(tau - t) + tau; % information

Ibox_dot = heaviside(tau - t);  % rate of information

% exponential filter

Fexp = -exp(lambda*t) * lambda; 

Iexp =  -2 * tau * (exp(lambda * t) - 1).^2 ./ (exp(2*lambda * t) - 1);

for i=1:length(t)
x = (exp(lambda * t(i)) - 1) / (exp(2*lambda * t(i))-1) * exp(lambda * t(i));
Iexp_dot(i) = 4 * x * (1 - x);
end

%% strong FF network filters

% strong FF flat filter

N = 10; % number of FF modes

for i=1:N
M(i,:) = (-lambda * t).^(N-i) / factorial(N-i);  % coefficients in filter equation
gam(i,:) = gamma(N-i+1) * gammainc( -lambda * t,  N-i+1, 'lower') /  factorial(N-i); % signal component for optimal readout weights

for j=1:N
   S(i,j,:) = gamma(2 * N-i-j+1) * gammainc(-2 * lambda * t, 2 * N-i-j+1, 'lower') / (2^(2*N - i - j + 1) * factorial(N-i) * factorial(N-j));  % covariance component for optimal readout weights
end
end

for i=1:length(t)
    wopt(:,i) = inv(S(:,:,i)) * gam(:,i);  % optimal readout weights (dynamic)
    IFF_opt(i) = gam(:,i)' * inv(S(:,:,i)) * gam(:,i);  % optimal information readout
    IFF_flat(i) =  sum(gam(:,i))^2 / sum(sum(S(:,:,i))); % flat decoder readout
end
for i=1:length(t)
        IFF_opt_longrun(i) = (wopt(:,end)' * gam(:,i))^2 / (wopt(:,end)' * S(:,:,i) * wopt(:,end));  % longrun optimal decoder
end

F_FFflat = exp(lambda * t) .* sum(M,1);  % flat filter
F_FF_opt_longrun =   exp(lambda * t) .* sum((wopt(:,end)).*M,1); % long run optimal filter
F_FF_opt =   exp(lambda * t) .* sum((wopt).*M,1);  % dynamically optimal filter

%% Now do a mixed FF/NonNormal Network
tau_a = 1;

qmax = 3;

Sigma_u = [1,0.0;0.0,1];

anglemax = 101;

tau_ratios = [0.1, 0.5, 0.75, 0.9];

thetadiff = [5, 10];

for k=1:2
for q=1:4

tau_b = tau_a * tau_ratios(q);
       
theta_u = pi/2;
theta_a = theta_u;

theta_b = theta_u + thetadiff(k) * pi/180;

va = [cos(theta_a); sin(theta_a)];
vb = [cos(theta_b); sin(theta_b)];
uprime = [cos(theta_u); sin(theta_u)];

normchoice = 'Vleftnorm';

if strcmp(normchoice, 'Vleftnorm')

    Vleft = [va'; vb'];
    Vright = inv(Vleft);

end

A = Vright * diag(-1./[tau_a, tau_b]) * Vleft;

for i=1:length(t)

Gamma(:,:,i) = Vleft * Sigma_u * Vleft' ./ [2 * tau_a , tau_a + tau_b; tau_a + tau_b, 2 * tau_b] .* ([(exp(-t(i)*(1/tau_a + 1/tau_a))), (exp(-t(i)*(1/tau_a + 1/tau_b))); (exp(-t(i)*(1/tau_a + 1/tau_b))), (exp(-t(i)*(1/tau_b + 1/tau_b)))]-ones(2)); 

IF(i,q,k) = -uprime' * (expm(A' * t(i)) - eye(2)) * Vleft' * inv(Gamma(:,:,i)) * Vleft * (expm(A * t(i)) - eye(2)) * uprime;

end

end
end

%% Plot results

% subplot 1: temporal filter

subplot(3,3,1)
hold on
plot(t,Fbox, 'linewidth', 3)
plot(t,Fexp, 'linewidth', 3)
set(gca, 'fontsize', 18)
xlabel('Time')
ylabel('Filter')
legend('Box Filter (width=\tau)', 'Exponential Filter (width=\tau)')
title('Temporal Filter')
axis([0,10,0,1])


subplot(3,3,2)
hold on
plot(t,F_FFflat, 'linewidth', 3, 'color', 'b')
plot(t,F_FF_opt_longrun, 'linewidth', 3, 'color', 'r')
plot(t,F_FF_opt, 'linewidth', 3, 'color', 'g')
set(gca, 'fontsize', 18)
xlabel('Time')
ylabel('Filter')
legend('Fixed Readout (Flat)', 'Fixed Readout (Long-Run Optimum)', 'Time-Varying Readout (Optimum)')
title('Response Projection (10 Strong FF Modes Network)')



% subplot 2: information vs time

subplot(3,3,4)
hold on
plot(t,Ibox, 'linewidth', 3)
plot(t,Iexp, 'linewidth', 3)
set(gca, 'fontsize', 18)
xlabel('Time')
ylabel('Information')
legend('Box Filter (width=\tau)', 'Exponential Filter (width=\tau)')
title('Output Information')
axis([0,10,0,2])

subplot(3,3,5)
hold on
plot(t,IFF_flat, 'linewidth', 3, 'color', 'b')
plot(t,IFF_opt_longrun, 'linewidth', 3, 'color', 'r')
plot(t,IFF_opt, 'linewidth', 3, 'color', 'g')
set(gca, 'fontsize', 18)
xlabel('Time')
ylabel('Information')
legend('Fixed Readout (Flat)', 'Fixed Readout (Long-Run Optimum)', 'Time-Varying Readout (Optimum)')
title('Network Response Information')

subplot(3,3,6)
clrs = {'k', [100,100,100]/255, [200,200,200]/255};
hold on
for i=1:3
plot(t,IF(:,i,1), 'linewidth', 3, 'color', clrs{i})
end
% for i=1:3
% plot(t,IF(:,i,2), 'linewidth', 3, 'linestyle', '--', 'color', clrs{i})
% end
set(gca, 'fontsize', 18)
xlabel('Time')
ylabel('Information')
legend('\tau_1 = 1, \tau_2 = 0.1', '\tau_1 = 1, \tau_2 = 0.5', '\tau_1 = 1, \tau_2 = 0.75')
title('Network Response Information')
axis([0,10,0,3.5])

% subplot 3: rate of information vs time

subplot(3,3,7)
hold on
plot(t,Ibox_dot, 'linewidth', 3)
plot(t,Iexp_dot, 'linewidth', 3)
set(gca, 'fontsize', 18)
xlabel('Time')
ylabel('Rate of Information')
legend('Box Filter (width=\tau)', 'Exponential Filter (width=\tau)')
axis([0,10,0,1])

subplot(3,3,8)
hold on
plot(t(1:(end-1)),diff(IFF_flat) ./ diff(t), 'linewidth', 3, 'color', 'b')
plot(t(1:(end-1)),diff(IFF_opt_longrun) ./ diff(t), 'linewidth', 3, 'color', 'r')
plot(t(1:(end-1)),diff(IFF_opt) ./ diff(t), 'linewidth', 3, 'color', 'g')
set(gca, 'fontsize', 18)
xlabel('Time')
ylabel('Rate of Information')
legend('Fixed Readout (Flat)', 'Fixed Readout (Long-Run Optimum)', 'Time-Varying Readout (Optimum)')

subplot(3,3,9)
hold on
for i=1:3
plot(t(1:(end-1)),diff(IF(:,i,1))' ./ diff(t), 'linewidth', 3, 'color', clrs{i})
end
% for i=1:3
% plot(t(1:(end-1)),diff(IF(:,i,2)) ./ diff(t), 'linewidth', 3, 'linestyle', '--', 'color', clrs{i})
% end
set(gca, 'fontsize', 18)
xlabel('Time')
ylabel('Rate of Information')
legend('\tau_1 = 1, \tau_2 = 0.1', '\tau_1 = 1, \tau_2 = 0.5', '\tau_1 = 1, \tau_2 = 0.75')
axis([0,10,0,1])