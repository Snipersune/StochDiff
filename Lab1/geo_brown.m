%% ====== Task 1 ====== %%
clear
close all

rng(42);

mu = 2; sigma = 1; X_0 = 1;

T = 10;
dt = 2^(-8);
t = dt:dt:T;
N = cast(floor(T/dt), "int16");

N_discs = 3;

sqrt_dt = sqrt(dt);
rands = randn(1,N);

dW = sqrt_dt*rands;
W = cumsum(dW);
X = exp((mu-0.5*sigma)*t + sigma*W);

figure();
plot(0:dt:T, [X_0,X], "DisplayName", sprintf("N=%d, dt=%0.4f", N, dt))
hold on
for i=2:N_discs
    dW = coarsen(dW, 2^2);
    N = size(dW,2);
    dt = T/N;
    t = dt:dt:T;
    W = cumsum(dW);
    X = exp((mu-0.5*sigma)*t + sigma*W);
    plot(0:dt:T, [X_0,X], "DisplayName", sprintf("N=%d, dt=%0.4f", N, dt))
    hold on
end
set(gca, "YScale", "log")
legend(Location="northwest")
xlabel('t', FontSize=14)
ylabel('(t)', FontSize=14)



%% ====== Task 2 ====== %%
clear
close all

rng(42)

mu = 2; sigma = 1; X_0 = 1;

T = 10;
dt = 2^(-9);
t = dt:dt:T;
N = cast(floor(T/dt), "int16");

sqrt_dt = sqrt(dt);

N_procs = [100, 1000, 10000];

means_Xt = zeros(size(N_procs,2),N);
vars_Xt = zeros(size(N_procs,2),N);

for i=1:size(N_procs,2)
    N_p = N_procs(i);
    rands = randn(N_p,N);
    
    dW = sqrt_dt*rands;
    W = cumsum(dW,2);
    
    X = exp((mu-0.5*sigma)*repmat(t,[N_p,1]) + sigma*W);

    means_Xt(i,:) = mean(X,1);
    vars_Xt(i,:) = var(X,1,1);
end


figure();
for i=1:size(N_procs,2)
    mean_Xt = means_Xt(i,:);
    plot(0:dt:T, [X_0,mean_Xt], "DisplayName", sprintf("Realizations=%d", N_procs(i)))
    hold on
end
legend()
set(gca, "YScale", "log")
xlabel('t', FontSize=14)
ylabel('E[X(t)]',FontSize=14)

figure();
for i=1:size(N_procs,2)
    var_Xt = vars_Xt(i,:);
    plot(0:dt:T, [0,var_Xt], "DisplayName", sprintf("Realizations=%d", N_procs(i)))
    hold on
end
legend()
set(gca, "YScale", "log")
xlabel('t', FontSize=14)
ylabel('V[X(t)]', FontSize=14)


%% ====== Task 3 ====== %%
clear
close all

rng(42)

mu = 2; sigma = 1; X_0 = 1;

T = 10;
dt = 2^(-9);
t = dt:dt:T;
N = cast(floor(T/dt), "int16");

sqrt_dt = sqrt(dt);

N_procs = 10000;
rands = randn(N_procs,N);
    
dW = sqrt_dt*rands;
W = cumsum(dW,2);
    
X = exp((mu-0.5*sigma)*repmat(t,[N_procs,1]) + sigma*W);
mean_Xt = mean(X,1);


figure();
plot(0:dt:T, [X_0,mean_Xt], "b-")
hold on
for i=1:5
    plot(0:dt:T, [X_0,X(i,:)], "r-.")
end    
legend(sprintf("E[X(t)] for %d realizations", N_procs), "5 individual paths", Location="northwest")
set(gca, "YScale", "log")
xlabel('t', FontSize=14)
ylabel('X(t)',FontSize=14)





