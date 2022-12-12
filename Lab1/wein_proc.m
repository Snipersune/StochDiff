%% ====== Task 1 ====== %%
clear
close all

rng(42);

T = 1; 
N = 512; 
dt = T/N;

rands = randn(1,N);
sqrt_dt = sqrt(dt);

dW = sqrt_dt.*rands;
W = cumsum(dW);

figure();
plot(0:dt:T, [0,W], "r-")
xlabel('t')
ylabel('W(t)')

%% ====== Task 2 ====== %%
clear
close all

rng(42)

T = 1; 
N = 256; 
dt = T/N;
N_discs = 4;

sqrt_dt = sqrt(dt);
rands = randn(1,N);

dW = sqrt_dt*rands;
W = cumsum(dW);

figure();
plot(0:dt:T, [0,W], "DisplayName", sprintf("N=%d, dt=%0.4f", N, dt))
hold on
for i=2:N_discs
    dW = coarsen(dW, 2);
    N = size(dW,2);
    dt = T/N;
    W = cumsum(dW);

    plot(0:dt:T, [0,W], "DisplayName", sprintf("N=%d, dt=%0.4f", N, dt))
    hold on
end
legend()
xlabel('t')
ylabel('W(t)')

%% ====== Task 3 ====== %%
clear
close all

rng(42)

T = 2;
dt = 2^(-7);
N = cast(floor(T/dt), "int16");

sqrt_dt = sqrt(dt);

N_procs = [100, 1000, 10000];

means_Wt = zeros(size(N_procs,2),N);
vars_Wt = zeros(size(N_procs,2),N);

for i=1:size(N_procs,2)
    N_p = N_procs(i);
    rands = randn(N_p,N);
    
    dW = sqrt_dt*rands;
    W = cumsum(dW,2);

    means_Wt(i,:) = mean(W,1);
    vars_Wt(i,:) = var(W,1,1);
end

figure();
for i=1:size(N_procs,2)
    mean_Wt = means_Wt(i,:);
    plot(0:dt:T, [0,mean_Wt], "DisplayName", sprintf("Realizations=%d", N_procs(i)))
    hold on
end
legend()
xlabel('t', FontSize=14)
ylabel('E[W(t)]',FontSize=14)

figure();
for i=1:size(N_procs,2)
    var_Wt = vars_Wt(i,:);
    plot(0:dt:T, [0,var_Wt], "DisplayName", sprintf("Realizations=%d", N_procs(i)))
    hold on
end
legend()
xlabel('t', FontSize=14)
ylabel('V[W(t)]', FontSize=14)




