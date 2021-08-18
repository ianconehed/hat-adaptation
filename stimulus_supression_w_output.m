clear all
close all

% load('u_it0.mat')

rng(9);
t_total = 2001;
dt = 1;
norm_noise = .2;
N = 1000; %number of neurons in LSM
tau_net = 10; %time constant for excitatory synaptic activation
g = 1.9; %LSM gain factor
p_c = .25;
r_0 = .2;
f_x = .0004;
f_y = .0002;
omega = 2*pi*(.004);
I = .0;
k_length = 2;


% J_ij = normrnd(0,g/sqrt(N),N,N); %LSM recurrent weights
J_ij = (g/sqrt(p_c*N))*full(sprandn(N,N,p_c));
O_ij = normrnd(0,g/sqrt(N),1,N);
ret_ij = normrnd(0,.3,N,1);
eig_W = eig(J_ij);
theta_i = rand(N,1)*2*pi;
PI_it_t1 = I*cos(omega*(1:t_total) + theta_i);
PI_it_x1 = zeros(N,t_total);
PI_it_x1(:,1:100) = I;
rate_it_record = zeros(N,t_total,k_length);

u_it0 = normrnd(0,1,N,1);

% figure
% plot(PI_it_t1(1,:))

for k = 1:k_length
    if k == 1
        PI_it_t = zeros(N,t_total); %input signal in t
%         PI_it_t = PI_it_t1;
        PI_it_x = zeros(N,t_total); %input signal in x
        PI_it_y = zeros(N,t_total); %input signal in y
    elseif k == 2
%         PI_it_t = PI_it_t1;
%         PI_it_t = abs(1*sin(2*pi*f_x*(1:t_total)));
%         PI_it_t(1:1000) = 1;
%         PI_it_t = zeros(1,t_total);
%         PI_it_x = PI_it_x1;
%         PI_it_x = (-2500:2500).^2; 
%         PI_it_x = PI_it_x./max(PI_it_x); 
%         PI_it_y = zeros(N,t_total); 
    end
    u_it = zeros(N,t_total);
    z = zeros(1,t_total);
    rate_it = zeros(N,t_total);
%     u_it(:,1) = normrnd(0,1,N,1);
    u_it(:,1) = u_it0;

    figure('Renderer', 'painters', 'Position', [200 300 1000 700])
    subplot(2,2,1)
    plot(PI_it_t(1,:))
    hold on
    plot(PI_it_x(1,:))
    hold on
    plot(PI_it_y(1,:))
    ylabel('inputs')
    xlabel('Time (ms)')

    for t = 2:t_total
        rate_it(:,t-1) = phi(u_it(:,t-1),r_0);
        sum_wji = J_ij*rate_it(:,t-1);
        sum_oji = O_ij*rate_it(:,t-1);
        sum_ret = ret_ij*z(t-1);
        position_input = PI_it_t(:,t) + PI_it_x(:,t) + PI_it_y(:,t);
        z_input = PI_it_t1(1,t);
        del_ui = (randn(N,1)*norm_noise + position_input -u_it(:,t-1) + sum_wji + sum_ret)*(dt/tau_net);
        del_z = (-z(t-1) + sum_oji + z_input)*(dt/tau_net);
        u_it(:,t) = u_it(:,t-1) + del_ui;
        z(t) = z(t-1) + del_z;
    end
    
    rate_it_record(:,:,k) = rate_it;
    
    a_corr = zeros(1,t_total);
    for tau = 1:t_total
        a_corr_tau = rate_it(:,1:(t_total-tau)).*rate_it(:,tau+1:t_total);
        a_corr_sum = sum(a_corr_tau,1)/N;
        a_corr(tau) = sum(a_corr_sum)/(t_total-tau);
    end
        
    subplot(2,2,2)
    imagesc(rate_it)
    ylabel('RNN activity')
    xlabel('Time (ms)')
    
    subplot(2,2,3)
    plot(a_corr)
    ylabel('Autocorrelation')
    xlabel('Time (ms)')
    
    subplot(2,2,4)
    plot(z)
    ylabel('output z')
    xlabel('Time (ms)')
end

corr_12 = zeros(1,t_total);
for t = 1:t_total
    mean_t = mean(rate_it(:,t));
    corr_12(t) = corr(rate_it_record(:,t,1),rate_it_record(:,t,2));
end

figure
plot(corr_12)
ylabel('Correlation betweeen state 1 and state 2 at time t')
xlabel('Time (ms)')


% figure
% plot(eig_W,'o')
% hold on
% plot(cos(0:.01:2*pi),sin(0:.01:2*pi));
% xlabel('Real component')
% ylabel('imaginary component')
% title('Eigenvalues')
% xlim([-1.5 1.5])
% ylim([-1.5 1.5])

function y = phi(x,r_0)
y = r_0*tanh(x/r_0).*(x<=0) + (2-r_0)*tanh(x/(2-r_0)).*(x>0);
end
