clearvars -except eye_old
close all

rng(8);
t_total = 5001;
dt = 1;
norm_noise = 0;
N = 1000; %number of neurons in LSM
tau_net = 50; %time constant for excitatory synaptic activation
tau_v = 10;
g = 1.7; %LSM gain factor
f_t = .0004;
p_c = .1;


W_ji = (g/sqrt(p_c*N))*full(sprandn(N,N,p_c)); %LSM recurrent weights
% W_ji = (W_ji + .1*W_ji.');
eig_W = eig(W_ji);
u_it0 = normrnd(0,2,N,1);


[W_in_t, W_in_x, W_in_y] = deal(normrnd(0,1,N,1)); %input weights

for k = 1:2
    PI_it_t = zeros(1,t_total); %input signal in t
    PI_it_x = zeros(1,t_total); %input signal in x
    PI_it_y = zeros(1,t_total); %input signal in y
    
    if k == 1
%         PI_it_t(1:1000) = 1;
    elseif k == 2
        PI_it_t = abs(1*sin(2*pi*f_t*(1:t_total)));
    end

    u_it = zeros(N,t_total);
    rate_it = zeros(N,t_total);
    u_it(:,1) = u_it0;


    for t = 2:t_total
        rate_it(:,t-1) = tanh(u_it(:,t-1));
        sum_wji = W_ji*rate_it(:,t-1);
        position_input = W_in_t*PI_it_t(t-1) + W_in_x*PI_it_x(t-1) + W_in_y*PI_it_y(t-1);
        del_ui = (randn(N,1)*norm_noise + position_input -u_it(:,t-1) + sum_wji)*(dt/tau_net);
        u_it(:,t) = u_it(:,t-1) + del_ui;
    end

    [max_am1,max_in1] = max(rate_it,[],2);
    [be1,eye1] = sort(max_in1);
    if k == 1
        eye_old = eye1;
        figure('Renderer', 'painters', 'Position', [200 300 1000 700])
        subplot(3,2,1)
        imagesc(rate_it(eye1,:))
        title('Network activity')
        
        subplot(3,2,3)
        imagesc(rate_it(eye1,:))
        title('Network activity')
    elseif k == 2
        subplot(3,2,2)
        imagesc(rate_it(eye_old,:))
        title('Network activity w/ stimulus, same sorting')
        subplot(3,2,4)
        imagesc(rate_it(eye1,:))
        title('Network activity w/ stimulus, new sorting')
    end

end


subplot(3,2,5)
imagesc(W_ji)
title('Weights')
subplot(3,2,6)
plot(eig_W,'o')
hold on
plot(cos(0:.01:2*pi),sin(0:.01:2*pi));
xlabel('Real component')
ylabel('imaginary component')
title('Eigenvalues')
xlim([-1.5 1.5])
ylim([-1.5 1.5])

