clearvars -except V_circ
close all

N = 50;
del_x = 1.87;
del_t = 6.1;
dt = .001;
v = del_x/del_t;
t = 0:dt:del_t - dt*del_t;
x = v*t;
sigma = .9/(3*sqrt(2));
alpha = 7.5;

[R_i] = deal(zeros(N,length(x)));
[LTP,LTD,LTP_W,LTD_W] = deal(zeros(N,length(t)));
[IP,ID] = deal(zeros(1,N));

tau_p = .2;
tau_d = 1.5;
eta_p = .2;
eta_d = 200;
% eta_W = .0002;
eta_W = .0002;
max_p = 2.2;
max_d = 2;
tau_P = .4;

p_mag = 1;

num_laps = 10;
num_tests = 1;

tr_loc = round((N)/2);

for i = 1:N
    R_i(i,:) = circshift(exp(-((x-(del_x/2))).^2/(2*(sigma^2))),i*round((length(x)/N)));
end

R_i = [R_i(N/2 + 1:N,:);R_i(1:N/2,:)];
% rand_mix = randperm(N);
% R_i = R_i(rand_mix,:);
% W = W(rand_mix);

V_diff_fin = zeros(num_tests,8003);
V_init_fin = zeros(num_tests,8003);

plastic = zeros(8003,1000);



for test = 1:num_tests
    W_center = randperm(N,1);
%     W_center = 4;
%     W = 0 + .17*exp(-(((1:N) - W_center).^2)/65);
    W = 0 + 0.0000000*.15*exp(-(((1:N) - W_center).^2)/15);
    
%     p_loc = randperm(N,1)*(del_x/N);
    p_loc = (N/2)*(del_x/N);
%     p_loc = 25*(del_x/N);
    p_time = p_loc/v;
    p_t_ind = round(p_time/dt);
    
    V_t = zeros(num_laps,length(t));
    V_t(1,:) = alpha*(W*R_i);
    
    P = zeros(1,length(t));
    
    
    
    
    if test == num_tests
        figure('Position',[0 0 2000 2000])
        subplot(6,3,1)
        imagesc(W)
        title('Weights before lap 1')
        colorbar
        subplot(6,3,2)
        plot(t,V_t,'k')
        hold on
        plot(t,10*LTP(tr_loc,:),'r--')
        hold on
        plot(t,5*LTD(tr_loc,:),'b--')
        hold on
        plot(t,5*P,'g')
        title('CA1 ramp depolarization')
        ylim([0 20])
        xlim([0 max(t)])
        xlabel('time(seconds)')
        ylabel('ramp amplitude(mV)')
        hold off
    end
    
    for lap = 2:num_laps
        dW = 0;
        for i = 2:length(t) %set to 1:length(t) to get looping
            if i == 1
                t_minus = length(t);
            else
                t_minus = i-1;
            end
            t_val = t(t_minus);
            x_val = v*t_val;
            x_minus = find(x == x_val);
            dLTP = (-LTP(:,t_minus) + eta_p*R_i(:,x_minus).*(max_p-LTP(:,t_minus)))*(dt/tau_p);
            dLTD = (-LTD(:,t_minus) + eta_d*R_i(:,x_minus).*(max_d-LTD(:,t_minus)))*(dt/tau_d);
            LTP(:,i) = LTP(:,t_minus) + dLTP;
            LTD(:,i) = LTD(:,t_minus) + dLTD;
            LTP_W(:,i) = LTP(:,i).*(1-W');
            LTD_W(:,i) = LTD(:,i).*(W');

            if i == round(p_time/dt)
                P(i) = p_mag;
            else
                P(i) = P(t_minus) - P(t_minus)*(dt/tau_P);
            end

            dW = dW + eta_W*P(i)*(LTP_W(:,i)-LTD_W(:,i));
        end
        if lap == 3 && test == num_tests
            IP = LTP*P';
            ID = LTD*P';
            W_fixed = (IP./(IP + ID))';
        end
        W = W + dW';%.*(1-W);
        W(W<0) = 0;
        
        V_t(lap,:) = alpha*(W*R_i);
        
        
        if test == (num_tests) && lap > num_laps-4
            lap2 = lap - (num_laps-5);
            subplot(6,3,((lap2-1)*3) + 1)
            imagesc(W)
            title(['Weights at lap ',num2str(lap)])
            colorbar
            subplot(6,3,((lap2-1)*3)+2)
            plot(t,V_t(lap,:),'k')
            hold on
            plot(t,10*LTP(tr_loc,:),'r--')
            hold on
            plot(t,5*LTD(tr_loc,:),'b--')
            hold on
            plot(t,5*P,'g')
            title('CA1 ramp depolarization')
            ylim([0 20])
            xlim([0 max(t)])
            xlabel('time(seconds)')
            ylabel('ramp amplitude(mV)')
            hold off
        end
    end
    
    p_t_ind2 = length(t) + round(p_time/dt);
    V_init = [V_t(1,:) V_t(1,:) V_t(1,:)];
    V_final = [V_t(num_laps-2,:) V_t(num_laps-1,:) V_t(num_laps,:)];
    
    V_diff = V_final-V_init;
    V_diff_fin(test,:) = V_diff(p_t_ind2-length(0:dt:4):p_t_ind2+length(0:dt:4));
    V_init_fin(test,:) = V_init(p_t_ind2-length(0:dt:4):p_t_ind2+length(0:dt:4));
    
    for k = 1:8003
        plastic(k,round((V_init_fin(test,k)+1)*100)) = V_diff_fin(test,k);
    end
end


ax(1) = subplot(2,3,3);
for i = 1:num_tests
    cline((1:length(V_diff_fin(i,:)))*dt-4,V_diff_fin(i,:),[],V_init_fin(i,:))
    hold on
end
hold off
ylabel('Change in ramp amplitude (mV)')
xlabel('Time relative to plateau onset (s)')
xlim([-4 4])
c = colorbar;
c.Label.String = 'Initial ramp amplitude (mV)';
colormap(ax(1),jet)


[x_int,y_int,z_int] = find(plastic);
x_int = (x_int-4001)*dt;
y_int = (y_int/100) - 1;
[xq,yq] = meshgrid(-4:.01:4,0:.01:10);
vq = griddata(x_int,y_int,z_int,xq,yq,'cubic');
vq = vq- 1.5;

iblur = imgaussfilt(vq,6);

ax(2) = subplot(2,3,6);
imagesc([-4 4], [0 10], vq)
set(gca,'YDir','normal')
ylabel('Initial ramp amplitude (mV)')
xlabel('Time relative to plateau onset (s)')
d = colorbar;
d.Label.String = 'Change in ramp amplitude (mV)';

colormap(ax(2),brewermap([],'*RdBu'))

subplot(6,3,16)
imagesc(W_fixed)
title('fixed point weights')
colorbar


% figure
% plot((1:N)*(del_t/N) - (del_t/2),W_fixed(end:-1:1))
% title("fixed point linear")
% xlabel("time relative to plateau onset")
% ylabel("synaptic strength")
% 
% 
V_new = alpha*W_fixed*R_i;
% state_fixed = 1 - V_t./V_new;
% counter = num_laps - sum(state_fixed < .1,1);
% % counter = counter(1:round(length(counter)/50):end);
% convergence_t = 1./(IP + ID);
% conv_t_ex = interp1(1:N,convergence_t,1:N/(length(counter)):N);
% figure
% yyaxis left
% plot(conv_t_ex)
% yyaxis right
% plot(counter)
% 
% counter_v_t = counter(1:length(conv_t_ex))./conv_t_ex;
% counter_v_t_dist = counter_v_t./(abs(V_t(1,length(conv_t_ex):-1:1)-V_new(length(conv_t_ex):-1:1)));
% figure
% scatter(counter(1:length(conv_t_ex)),conv_t_ex)
% 
% 
% figure
% imagesc((1-V_t./V_new)*100)
% set(gca,'ColorScale','log')
% xlabel('timesteps')
% ylabel('trials')
% title('Percent away from fixed point')
% 
figure
colororder(jet(5))
plot(t - max(t)/2,V_t(1,end:-1:1))
hold on
plot(t - max(t)/2,V_t(2,end:-1:1),'--')
hold on
plot(t - max(t)/2,V_t(4,end:-1:1),'--')
hold on
plot(t - max(t)/2,V_t(6,end:-1:1),'--')
% hold on
% plot(W_fixed(end:-1:1),'k-')
title("fixed point linear")
xlabel("time relative to plateau onset")
ylabel("voltage")
% xlim([round(-del_t/2) round(del_t/2)])
% ylim([0 14])

% figure
% subplot(1,3,1)
% imagesc(R_i)
% title('CA3 inputs')
% subplot(1,3,2)
% imagesc(W)
% title('Weights')
% subplot(1,3,3)
% imagesc(V_i)
% title('CA1 ramp depolarization')

% figure
% imagesc(V_i)
% title('CA1 ramp depolarization')
W_fixed_old = W_fixed;