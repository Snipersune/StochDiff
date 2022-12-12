clear all 
close all

DIM = 2;
rng(42);

T = 1; 
N = 1024; 
dt = T/N;

n_samples = 1000;

rands = randn(n_samples,N);
sqrt_dt = sqrt(dt);

dW = sqrt_dt.*rands;
W = cumsum(dW,2);

ito_I1 = sum(W(:,1:end-1)*dt, DIM);
strat_I1 = sum((0.5*(W+[zeros(n_samples,1),W(:,1:end-1)]) + 0.5*sqrt_dt*randn(n_samples,N))*dt, DIM);

ito_I2 = sum([zeros(n_samples,1),W(:,1:end-1)].*dW, DIM);
strat_I2 = sum((0.5*(W+[zeros(n_samples,1),W(:,1:end-1)]) + 0.5*sqrt_dt*randn(n_samples,N)).*dW, DIM);

ito_I2_true = 0.5*(W(:,end).^2-T);
strat_I2_true = 0.5*W(:,end).^2;

ito_I2_err = abs(ito_I2 - ito_I2_true);
strat_I2_err = abs(strat_I2 - strat_I2_true);



figure()
histogram(ito_I2_err)
ylabel("Error")
legend("Ito")

figure()
histogram(strat_I2_err)
label("Error")
legend("Strat")





