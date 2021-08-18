clear all
close all

% load('u_it0.mat')

% rng(5);
t_total = 5001;
dt = 1;
norm_noise = .1;
N = 1000; %number of neurons in LSM
M = 100; %output neurons
tau_net = 50; %time constant for excitatory synaptic activation
tau_v = 50;
g = 1.7; %LSM gain factor
p_c = .1;
f_x = .0004;
f_y = .0002;


W_ji = (g/sqrt(p_c*N))*full(sprandn(N,N,p_c)); %LSM recurrent weights
O_jk = (g/sqrt(p_c*N))*full(sprandn(N,M,p_c)); %output weights
[W_in_t, W_in_x, W_in_y] = deal(10*normrnd(0,g/sqrt(N),N,1)); %input weights
eig_W = eig(W_ji);
u_it0 = normrnd(0,1,N,1);
% PI_it(:,2501:3500) = repmat(-.001*(1:1000) + 1,N,1);
% PI_it(:,3501:7001) = repmat(PI_it(:,3500),1,3501);
% PI_it(:,7002:15001) = PI_it(:,1:8000);
% PI_it(:,15001:35001) = 0;

for k = 1:2
    if k == 1
        PI_it_t = zeros(1,t_total); %input signal in t
        PI_it_x = zeros(1,t_total); %input signal in x
        PI_it_y = zeros(N,t_total); %input signal in y
        PI_it_t = abs(2*sin(2*pi*f_x*(1:t_total)));
        PI_it_y(:,1:1000) = 2*ones(N,1000).*normrnd(0,1,N,1);
%         PI_it_y(:,1:1000) = .2;
    elseif k == 2
        PI_it_y(:,1:1000) = 2*ones(N,1000).*normrnd(0,1,N,1);
%         PI_it_t = abs(1*sin(2*pi*f_x*(1:t_total)));
%         PI_it_t(1:1000) = 1;
%         PI_it_t = zeros(1,t_total);
%         PI_it_x = zeros(1,t_total);
%         PI_it_x = (-2500:2500).^2; 
%         PI_it_x = PI_it_x./max(PI_it_x); 
%         PI_it_y = zeros(1,t_total); 
    elseif k == 3
%         PI_it_t = abs(1*sin(2*pi*f_x*(1:t_total)));
%         PI_it_t(1:1000) = 1;
%         PI_it_x = (-2500:2500).^2; 
%         PI_it_x = PI_it_x./max(PI_it_x); 
%         PI_it_y = zeros(1,t_total);
    elseif k == 4
%         PI_it_t = abs(1*sin(2*pi*f_x*(1:t_total)));
%         PI_it_t(1:1000) = 1;
%         PI_it_x = zeros(1,t_total);
%         PI_it_x = (-2500:2500).^2; 
%         PI_it_x = PI_it_x./max(PI_it_x); 
%         PI_it_y = 0*sin(2*pi*f_y*(1:t_total)); %input position signal in y
%         PI_it_y(:,1:5000) = .0002*(1:5000);
    end
    u_it = zeros(N,t_total);
    rate_it = zeros(N,t_total);
%     u_it(:,1) = normrnd(0,2,N,1);
    u_it(:,1) = u_it0;
    v_kt = zeros(M,t_total);
    rate_kt = zeros(M,t_total);

    figure('Renderer', 'painters', 'Position', [200 300 1000 700])
    subplot(2,2,1)
    plot(PI_it_t)
    hold on
    plot(PI_it_x)
    hold on
    plot(PI_it_y(1,:))
    ylabel('inputs')
    xlabel('Time (ms)')

    for t = 2:t_total
        rate_it(:,t-1) = tanh(u_it(:,t-1));
        rate_kt(:,t-1) = tanh(v_kt(:,t-1)-1) + 1;
        sum_wji = W_ji*rate_it(:,t-1);
        sum_output = O_jk'*rate_it(:,t-1);
        position_input = W_in_t*PI_it_t(t-1) + W_in_x*PI_it_x(t-1) + PI_it_y(t-1);
        del_ui = (randn(N,1)*norm_noise + position_input -u_it(:,t-1) + sum_wji)*(dt/tau_net);
        u_it(:,t) = u_it(:,t-1) + del_ui;

        del_vk = (sum_output -v_kt(:,t-1))*(dt/tau_v);
        v_kt(:,t) = v_kt(:,t-1) + del_vk;
    end

    [max_am,max_in] = max(rate_kt,[],2);
    [be,eye] = sort(max_in);
    if k ==1
        eye_orig = eye;
    end
    v_kt_sorted_orig = rate_kt(eye_orig,:);
    v_kt_sorted = rate_kt(eye,:);
    
    [max_am1,max_in1] = max(rate_it,[],2);
    [be1,eye1] = sort(max_in1);
    

    subplot(2,2,2)
    imagesc(rate_it(eye1,:))
    ylabel('RNN activity')
    xlabel('Time (ms)')

    subplot(2,2,3)
    imagesc(v_kt_sorted)
    ylabel('output activity w/ new sorting')
    xlabel('Time (ms)')
    
    subplot(2,2,4)
    imagesc(v_kt_sorted_orig)
    ylabel('output activity w/ original sorting')
    xlabel('Time (ms)')
end

figure
plot(eig_W,'o')
hold on
plot(cos(0:.01:2*pi),sin(0:.01:2*pi));
xlabel('Real component')
ylabel('imaginary component')
title('Eigenvalues')
xlim([-1.5 1.5])
ylim([-1.5 1.5])