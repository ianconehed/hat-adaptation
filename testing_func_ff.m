clear all
close all

rng(9);
t_total = 10001;
dt = 1;
N = 100; %number of neurons in LSM
tau_net = 50; %time constant for excitatory synaptic activation

v = ones(N-1,1);
T = diag(v,1);
U_init = rand(N,N)-.5;
[U,R] = mgsog(U_init);
W_ji = U*T*U^(-1);
% W_ji = normrnd(0,g/sqrt(N),N,N); %LSM recurrent weights
% W_ji = (W_ji + .2*W_ji.')/2;
eig_W = eig(W_ji);

u_it = zeros(N,t_total);
% u_it(:,1) = U(1,:);
u_it(:,1) = normrnd(0,1,N,1);


for t = 2:t_total
    sum_wji = W_ji*u_it(:,t-1);
    del_ui = (-u_it(:,t-1) + sum_wji)*(dt/tau_net);
    u_it(:,t) = u_it(:,t-1) + del_ui;
end

proj = flip((U.'*u_it),1);


figure('Renderer', 'painters', 'Position', [200 300 1000 700])
subplot(2,1,1)
imagesc(u_it)
title('Network activity')
subplot(2,2,3)
imagesc(proj)
title('projection')
subplot(2,2,4)
plot(eig_W,'o')
hold on
plot(cos(0:.01:2*pi),sin(0:.01:2*pi));
xlabel('Real component')
ylabel('imaginary component')
title('Eigenvalues')
xlim([-1.5 1.5])
ylim([-1.5 1.5])