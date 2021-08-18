clear all
% close all

% rng(2);
t_total = 10001;
dt = 1;
N = 200; %number of neurons in LSM
tau_net = 50; %time constant for excitatory synaptic activation
norm_noise = 0.1;
f_t = .0004;
g = 1; %LSM gain factor
p_c = .1;
% (g/sqrt(p_c*N))*full(sprandn(N,N,p_c))

v = 1*ones(N-1,1);
T = 0*diag(v,1) + 0*(g/sqrt(p_c*N))*full(sprandn(N,N,p_c));
% T(N/2,(N/2) + 1) = 0;
% rand_b = randperm(length(T));
% T2 = T(rand_b,rand_b);
% T1 = (T1 + T1.')/2;
% T1 = tril(T1);
U_init = rand(N,N)-.5;
[U,R] = mgsog(U_init);
rand_a = randperm(length(T));
T1 = T(rand_a,rand_a);

% W_ji = U*T*U^(-1);
% W_ji = U*T*U^(-1) + .0*normrnd(0,g/sqrt(N),N,N);
W_ji = 0*U*T*U^(-1) +  1.8*(g/sqrt(p_c*N))*full(sprandn(N,N,p_c));
% W_ji = (U*T1*U^(-1) + U*T2*U^(-1))/2;
% W_ji = normrnd(0,g/sqrt(N),N,N); %LSM recurrent weights
% W_ji = T1;
[W_in_t, W_in_x, W_in_y] = deal(normrnd(0,1/sqrt(N),N,1));


[eigvect,eigval] = eig(W_ji);
gamma = (.05*rand(N,N) + .95)./abs(eigval);
gamma(isinf(gamma)|isnan(gamma)) = 0;
eigval_normed = eigval.*gamma;
W_ji2 = real(eigvect*eigval_normed*inv(eigvect));
W_ji3 = real(eigvect*eigval*inv(eigvect));

[eigvect2,eigval2] = eig(W_ji2);

u_it = zeros(N,t_total);
rate_it = zeros(N,t_total);
% u_it(:,1) = U(1,:);
% u_it(1,1:10) = 1;
u_it(:,1) = normrnd(0,1,N,1);

PI_it_t = zeros(1,t_total); %input signal in t
PI_it_x = zeros(1,t_total); %input signal in x
PI_it_y = zeros(1,t_total); %input signal in y

% PI_it_t(1:1000) = 1;
% PI_it_t = abs(1*sin(2*pi*f_t*(1:t_total)));
rng(9)
for k = 1:2
    if k == 2
%         PI_it_t(1:1000) = 1;
%         PI_it_t = 1*abs(1*sin(2*pi*f_t*(1:t_total)));
    end
    for t = 2:t_total
        rate_it(:,t-1) = tanh(u_it(:,t-1));
        sum_wji = W_ji*rate_it(:,t-1);
        position_input = W_in_t*PI_it_t(t-1) + W_in_x*PI_it_x(t-1) + W_in_y*PI_it_y(t-1);
        del_ui = (randn(N,1)*norm_noise + position_input -u_it(:,t-1) + sum_wji)*(dt/tau_net);
        u_it(:,t) = u_it(:,t-1) + del_ui;
    end

    proj = flip((U.'*rate_it),1);
    
    if k == 2
        eye_old = eye1;
    end
    
    [max_am1,max_in1] = max(rate_it,[],2);
    [be1,eye1] = sort(max_in1);
    
    figure('Renderer', 'painters', 'Position', [200 300 1000 700])
    subplot(2,2,1)
    plot(PI_it_t)
    hold on
    plot(PI_it_x)
    hold on
    plot(PI_it_y)
    title('Inputs')
    subplot(2,2,2)
    imagesc(rate_it(eye1,:))
    title('Network activity')
    subplot(2,2,3)
    imagesc(proj)
    title('projection')
    if k == 1
        subplot(2,2,4)
        plot(eigval,'o')
        hold on
        plot(cos(0:.01:2*pi),sin(0:.01:2*pi));
        xlabel('Real component')
        ylabel('imaginary component')
        title('Eigenvalues')
        xlim([-1.5 1.5])
        ylim([-1.5 1.5])
    elseif k ==2
        subplot(2,2,4)
        imagesc(rate_it(eye_old,:))
        title('Network activity old sorting')
    end
end
