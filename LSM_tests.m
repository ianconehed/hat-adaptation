clear all
close all

% load('u_it0.mat')

rng(9);
t_total = 5001;
dt = 1;
norm_noise = 0;
N = 1000; %number of neurons in LSM
M = 100; %output neurons
k_max = 40;
tau_net = 50; %time constant for excitatory synaptic activation
g = 1.02; %LSM gain factor
f_x = .0004;
f_y = .0002;

W_ji = normrnd(0,g/sqrt(N),N,N); %LSM recurrent weights
O_jk = normrnd(0,g/sqrt(N),N,M); %output weights
[W_in_t, W_in_x, W_in_y] = deal(normrnd(0,g/sqrt(N),N,1)); %input weights
w_eig = eig(W_ji);


u_it0 = normrnd(0,2,N,1);
u_it_all = zeros(N,t_total,k_max);




for k = 1:k_max
%     PI_it_t = abs(1*sin(2*pi*f_x*(1:t_total))); %input signal in t
    PI_it_t = zeros(1,t_total); %input signal in t
    PI_it_x = zeros(1,t_total); %input signal in x
    PI_it_y = zeros(1,t_total); %input signal in y
    
    u_it = zeros(N,t_total);
    rate_it = zeros(N,t_total);
    u_it(:,1) = normrnd(0,2,N,1);
    
    for t = 2:t_total
        sum_wji = W_ji*u_it(:,t-1);
        position_input = W_in_t*PI_it_t(t-1) + W_in_x*PI_it_x(t-1) + W_in_y*PI_it_y(t-1);
        del_ui = (randn(N,1)*norm_noise + position_input -u_it(:,t-1) + sum_wji)*(dt/tau_net);
        u_it(:,t) = u_it(:,t-1) + del_ui;
        rate_it(:,t) = tanh(u_it(:,t));
    end
    u_it_all(:,:,k) = u_it;
end

u_it_mean1 = mean(u_it_all(:,:,1:(k_max/2)),3);
u_it_mean2 = mean(u_it_all(:,:,((k_max/2) + 1):k_max),3);
[max_am1,max_in1] = max(u_it_mean1,[],2);
[be1,eye1] = sort(max_in1);
u_it_sorted1 = u_it_mean1(eye1,:);
u_it_sorted2 = u_it_mean2(eye1,:);


figure('Renderer', 'painters', 'Position', [200 300 1000 700])
subplot(2,2,1)
imagesc(u_it_all(eye1,:,1))
xlabel('Time(ms)')
ylabel('Neuron number')
title('First trial')


subplot(2,2,2)
imagesc(u_it_all(eye1,:,2))
xlabel('Time(ms)')
ylabel('Neuron number')
title('Second trial')

subplot(2,2,3)
imagesc(u_it_sorted1)
xlabel('Time(ms)')
ylabel('Neuron number')
title('Average of first k_max/2 trials')

subplot(2,2,4)
imagesc(u_it_sorted2)
xlabel('Time(ms)')
ylabel('Neuron number')
title('Average of second k_max/2 trials, using same sorting as first k_max/2 average')


figure
plot(w_eig,'o')
hold on
plot(cos(0:.01:2*pi),sin(0:.01:2*pi));
xlabel('Real component')
ylabel('imaginary component')
xlim([-1.5 1.5])
ylim([-1.5 1.5])