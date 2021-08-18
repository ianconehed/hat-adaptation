clear all
close all

% load('u_it0.mat')

rng(9);
t_total = 5001;
dt = 1;
norm_noise = 0;
N = 1000; %number of neurons in LSM
M = 100; %output neurons
tau_net = 50; %time constant for excitatory synaptic activation
tau_v = 50;
g = 1.02; %LSM gain factor
f_x = .0004;
f_y = .0002;


W_ji = normrnd(0,g/sqrt(N),N,N); %LSM recurrent weights
O_jk = normrnd(0,g/sqrt(N),N,M); %output weights
[W_in_t, W_in_x, W_in_y] = deal(normrnd(0,g/sqrt(N),N,1)); %input weights
% PI_it(:,2501:3500) = repmat(-.001*(1:1000) + 1,N,1);
% PI_it(:,3501:7001) = repmat(PI_it(:,3500),1,3501);
% PI_it(:,7002:15001) = PI_it(:,1:8000);
% PI_it(:,15001:35001) = 0;

for k = 1:4
    if k == 1
        PI_it_t = zeros(1,t_total); %input signal in t
        PI_it_x = zeros(1,t_total); %input signal in x
        PI_it_y = zeros(1,t_total); %input signal in y
    elseif k == 2
        PI_it_t = abs(1*sin(2*pi*f_x*(1:t_total)));
%         PI_it_t(1:1000) = 1;
%         PI_it_t = zeros(1,t_total);
        PI_it_x = zeros(1,t_total);
%         PI_it_x = (-2500:2500).^2; 
%         PI_it_x = PI_it_x./max(PI_it_x); 
        PI_it_y = zeros(1,t_total); 
    elseif k == 3
        PI_it_t = abs(1*sin(2*pi*f_x*(1:t_total)));
%         PI_it_t(1:1000) = 1;
        PI_it_x = (-2500:2500).^2; 
        PI_it_x = PI_it_x./max(PI_it_x); 
        PI_it_y = zeros(1,t_total);
    elseif k == 4
        PI_it_t = abs(1*sin(2*pi*f_x*(1:t_total)));
%         PI_it_t(1:1000) = 1;
%         PI_it_x = zeros(1,t_total);
%         PI_it_x = (-2500:2500).^2; 
        PI_it_x = PI_it_x./max(PI_it_x); 
        PI_it_y = 0*sin(2*pi*f_y*(1:t_total)); %input position signal in y
        PI_it_y(:,1:5000) = .0002*(1:5000);
    end
    u_it = zeros(N,t_total);
    u_it(:,1) = normrnd(0,2,N,1);
    v_kt = zeros(M,t_total);

    figure('Renderer', 'painters', 'Position', [200 300 1000 700])
    subplot(2,2,1)
    plot(PI_it_t)
    hold on
    plot(PI_it_x)
    hold on
    plot(PI_it_y)
    ylabel('inputs')
    xlabel('Time (ms)')

    for t = 2:t_total
        sum_wji = W_ji*u_it(:,t-1);
        sum_output = O_jk'*u_it(:,t-1);
        position_input = W_in_t*PI_it_t(t-1) + W_in_x*PI_it_x(t-1) + W_in_y*PI_it_y(t-1);
        del_ui = (randn(N,1)*norm_noise + position_input -u_it(:,t-1) + sum_wji)*(dt/tau_net);
        u_it(:,t) = u_it(:,t-1) + del_ui;

        del_vk = (sum_output -v_kt(:,t-1))*(dt/tau_v);
        v_kt(:,t) = v_kt(:,t-1) + del_vk;
    end

    [max_am,max_in] = max(v_kt,[],2);
    [be,eye] = sort(max_in);
    if k ==1
        eye_orig = eye;
    end
    v_kt_sorted_orig = v_kt(eye_orig,:);
    v_kt_sorted = v_kt(eye,:);
    
    [max_am1,max_in1] = max(u_it,[],2);
    [be1,eye1] = sort(max_in1);
    

    subplot(2,2,2)
    imagesc(u_it(eye1,:))
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