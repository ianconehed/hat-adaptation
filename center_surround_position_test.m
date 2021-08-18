clear all
close all

% load('u_it0.mat')

rng(40);
t_total = 2001;
dt = 1;
norm_noise = .0;
N = 30^2; %number of neurons in LSM
tau_net = 10; %time constant for excitatory synaptic activation
g = .2; %LSM gain factor
p_c = .25;
r_0 = .2;
f_x = .0004;
f_y = .0002;
omega = 2*pi*(.004);
I = .2;
k_length = 2;

k_adap = .5;
tau_adap = 400;


%% creation of weight matrix
alpha = 1.05;
sigma = .6; %weird behavior under one
alpha_neg = .05;
sigma_neg = 3;
offset = 5;
alpha_in = 3;
sigma_in = 2;
beta = -13.5;
ind = 1:N;
sz = [sqrt(N) sqrt(N)];
[row,col] = ind2sub(sz,ind);
dist = zeros(N,N);
for i = 1:N
    for j = 1:N
        row_dist = row(i) - row(j);
        if row_dist > sqrt(N)/2
            row_dist = row_dist - sqrt(N);
        end
        col_dist = col(i) - col(j);
        if col_dist > sqrt(N)/2
            col_dist = col_dist - sqrt(N);
        end
        dist(i,j) = sqrt(row_dist^2 + col_dist^2);
%         dist(i,j) = dist(i,j) - sqrt(N)*int(dist(i,j)/sqrt(N));
    end
end
% J_ij = alpha*exp(-(dist.^2)/(2*(sigma^2)));
J_ij = alpha*exp(-(dist.^2)/(2*(sigma^2)))...
    - alpha_neg*exp(-((dist+offset).^2)/(2*(sigma_neg^2)))...
    - alpha_neg*exp(-((dist-offset).^2)/(2*(sigma_neg^2)));
I_ij = alpha_in*exp(-(dist.^2)/(2*(sigma_in^2)));
% J_ij(J_ij<.0001 & J_ij>1*(10^-200)) = beta/N;

% J_ij = J_ij + (g/sqrt(N))*normrnd(0,1,N,N); %LSM recurrent weights
J_ij = J_ij + (g/sqrt(p_c*N))*full(sprandn(N,N,p_c)); %LSM recurrent weights
J_ij = J_ij - diag(diag(J_ij));
% J_ij = normrnd(0,g/sqrt(N),N,N); %LSM recurrent weights
% J_ij = (g/sqrt(p_c*N))*full(sprandn(N,N,p_c));
eig_W = eig(J_ij);
rate_it_record = zeros(N,t_total,k_length);

%% input
in_length = 800;
a = 1:sqrt(N);
b = 10*ones(1,sqrt(N));
indin = sub2ind(sz,1:sqrt(N),10*ones(1,sqrt(N)));
indin2 = repelem(indin,round(in_length/sqrt(N)));
indin_final = indin2(1:in_length);


u_it0 = normrnd(0,.1,N,1);
% rng(50);
for k = 1:k_length
    if k == 1
        PI_it_t = zeros(N,t_total); %input signal in t
%         PI_it_t = PI_it_t1;
        PI_it_x = zeros(N,t_total); %input signal in x
%         PI_it_x = ((-1000:1000).^2).*normrnd(0,1,N,1)/(1000^2);
        PI_it_x(:,1001:1800) = 1*ones(N,in_length).*I_ij(:,indin_final);
%         PI_it_x(:,1:450) = 1*ones(N,450).*I_ij(:,randi(N));
%         PI_it_x(:,1:500) = 1*ones(N,500).*I_ij(:,randi(N));
        PI_it_y = zeros(N,t_total); %input signal in y
%         PI_it_x(:,1101:1200) = 1*ones(N,100).*I_ij(:,randi(N));
        PI_it_x(:,1:200) = .1*ones(N,200).*normrnd(0,1,N,1);
    elseif k == 2
        PI_it_x(:,1:200) = .1*ones(N,200).*normrnd(0,1,N,1);
%         PI_it_t = PI_it_t1;
%         PI_it_t(:,1:100) = 1*ones(N,100).*normrnd(0,1,N,1);
%         PI_it_x(:,501:600) = PI_it_x(:,1101:1200);
%         PI_it_x(:,1:100) = 1*ones(N,100).*I_ij(:,randi(N));
%         PI_it_x(:,1101:1200) = 0*ones(N,100).*I_ij(:,200);
%         PI_it_x(:,1101:1200) = 1*ones(N,100).*I_ij(:,randi(N));
%         PI_it_x(:,1:2000) = .1*ones(N,2000).*normrnd(0,1,N,1);
%         PI_it_y = zeros(N,t_total); %input signal in y
%         PI_it_x = zeros(N,t_total); %input signal in x
%         PI_it_t = abs(1*sin(2*pi*f_x*(1:t_total)));
%         PI_it_t = zeros(1,t_total);
%         PI_it_x = PI_it_x1;
%         PI_it_x = (-2500:2500).^2; 
%         PI_it_x = PI_it_x./max(PI_it_x); 
%         PI_it_y = zeros(N,t_total); 
    end
    u_it = zeros(N,t_total);
    adap_it = zeros(N,t_total);
    rate_it = zeros(N,t_total);
    patt_time = zeros(sqrt(N),sqrt(N),t_total);
    patty = zeros(sqrt(N),sqrt(N),t_total);
%     u_it(:,1) = normrnd(0,.1,N,1);
%     u_it(1,1:100) = 10;
%     u_it(:,1) = u_it0;

    figure('Renderer', 'painters', 'Position', [200 300 1000 400])
    subplot(1,2,1)
    plot(mean(PI_it_t,1))
    hold on
    plot(mean(PI_it_x,1))
    hold on
    plot(mean(PI_it_y,1))
    ylabel('inputs')
    xlabel('Time (ms)')
    
    
    if k == 1
        v = VideoWriter('movie','Indexed AVI');
        v.FrameRate = 240;
        v.Colormap = parula;
        open(v);
    end

    for t = 2:t_total
        rate_it(:,t-1) = phi(u_it(:,t-1),r_0);
        sum_wji = J_ij*rate_it(:,t-1);
        position_input = PI_it_t(:,t) + PI_it_x(:,t) + PI_it_y(:,t);
        del_ui = (randn(N,1)*norm_noise + position_input -u_it(:,t-1) - adap_it(:,t-1) + sum_wji)*(dt/tau_net);
        del_adap = (-adap_it(:,t-1) + k_adap*u_it(:,t-1))*(dt/tau_adap);
        u_it(:,t) = u_it(:,t-1) + del_ui;
        adap_it(:,t) = adap_it(:,t-1) + del_adap;
        for i = 1:N
            patt_time(row(i),col(i),t) = rate_it(i,t-1);
            patty(row(i),col(i),t) = PI_it_x(i,t-1);
        end
        if k == 1
            im = mat2gray(patt_time(:,:,t), [-.1 .75]);
            [im2,cmap] = gray2ind(im,64);
            writeVideo(v,im2);  
        end
    end
    close(v)
    
    rate_it_record(:,:,k) = rate_it;
    
        
    subplot(1,2,2)
    imagesc(rate_it)
    ylabel('RNN activity')
    xlabel('Time (ms)')
    
%     subplot(2,2,3)
%     plot(a_corr)
%     ylabel('Autocorrelation')
%     xlabel('Time (ms)')
end

% close(v2);



corr_12 = zeros(1,t_total);
for t = 1:t_total
    mean_t = mean(rate_it(:,t));
    corr_12(t) = corr(rate_it_record(:,t,1),rate_it_record(:,t,2));
end

% figure
% plot(corr_12)
% ylabel('Correlation betweeen state 1 and state 2 at time t')
% xlabel('Time (ms)')

corr_22 = zeros(N,t_total);
for t = 1:t_total
    mean_t = mean(rate_it(:,t));
    corr_22(:,t) = rate_it_record(:,t,1).*rate_it_record(:,t,2);
end

figure
imagesc(corr_22)
ylabel('neuronxneuron Correlation betweeen state 1 and state 2 at time t')
xlabel('Time (ms)')

[max_am1,max_in1] = max(rate_it_record(:,1:1000,1),[],2);
[be1,eye1] = sort(max_in1);
rate_it_sorted = rate_it_record(eye1,:,1);
rate_it2_sorted = rate_it_record(eye1,:,2);

figure
subplot(1,2,1)
imagesc(rate_it_sorted)
subplot(1,2,2)
imagesc(rate_it2_sorted)
% figure
% plot(eig_W,'o')
% hold on
% plot(cos(0:.01:2*pi),sin(0:.01:2*pi));
% xlabel('Real component')
% ylabel('imaginary component')
% title('Eigenvalues')
% xlim([-1.5 1.5])
% ylim([-1.5 1.5])
% 
% figure
% histogram(randn(N,1)*norm_noise)
% hold on
% histogram(patty(:,:,2))


function y = phi(x,r_0)
y = r_0*tanh(x/r_0).*(x<=0) + (2-r_0)*tanh(x/(2-r_0)).*(x>0);
end
