clear all 
close all

DIM = 2;
rng(42);

T = 1; 
N = 256; 
dt = T/N;
ts = 0:dt:T;

n_samples = 10^6;

rands = randn(n_samples,N);
sqrt_dt = sqrt(dt);

dW = sqrt_dt.*rands;
W = cumsum(dW,2);

f1 = @(s, Ws) cos(s);
f2 = @(s, Ws) s*Ws;
f3 = @(s, Ws) s^3*sqrt(Ws);

I_f1 = sum(cos(ts(1:end-1)).*dW, DIM);
I_f2 = sum(ts(1:end-1).*[zeros(n_samples,1), W(:,1:end-1)].*dW, DIM);
I_f3 = sum(ts(1:end-1).^3.*sqrt(abs([zeros(n_samples,1), W(:,1:end-1)])).*dW, DIM);

mean_I_f1 = mean(I_f1)
mean_I_f2 = mean(I_f2)
mean_I_f3 = mean(I_f3)
%%
figure()
boxchart([I_f1, I_f2, I_f3])
hold on
yline(0, "r")

%%

figure()
histogram(I_f1)

figure()
histogram(I_f2)

figure()
histogram(I_f3)

