clear all
close all


a = .3;
k = .5;
mpa = 6;
spa = 5;
spa_n = 0;
alpha = .7;
sigma = .7;
alpha_neg = .1;
sigma_neg = 8;
offset = 7;
offest_zero = 20;

x = linspace(-50,50,10000);
y1 = mpa*cos(a*x).*(abs(x)<=pi/(2*a)) ...
    + spa*cos((x + ((k-1)*pi/(2*a)))/(k/a)).*(cos((x + ((k-1)*pi/(2*a)))/(k/a))>0).*(x> pi/(2*a))...
    + spa_n*cos((x + ((k-1)*pi/(2*a)))/(k/a)).*(cos((x + ((k-1)*pi/(2*a)))/(k/a))<0).*(x> pi/(2*a))...
    + spa*cos((x - ((k-1)*pi/(2*a)))/(k/a)).*(cos((x - ((k-1)*pi/(2*a)))/(k/a))>0).*(x< -pi/(2*a))...
    + spa_n*cos((x - ((k-1)*pi/(2*a)))/(k/a)).*(cos((x - ((k-1)*pi/(2*a)))/(k/a))<0).*(x< -pi/(2*a));

y1 = y1.*(abs(x)<=(2*pi)*(1 + 1/k));

y2 = alpha*exp(-(x.^2)/(2*(sigma^2)))...
    - alpha_neg*exp(-((x+offset).^2)/(2*(sigma_neg^2)))...
    - alpha_neg*exp(-((x-offset).^2)/(2*(sigma_neg^2)));


% f = @(x) (-1).*((x < -1)) + (x.^2+x).*((-1 <= x) & (x < 2)) + (x+2).*((2 <= x) & (x < 4)) + (x.^3+2*x.^2-3).*((4 <= x) & (x < 6)) + (2*exp(0.1*x-0.6)).*((6 <= x));
% x = linspace(-5, 10, 500);
% figure(1)
% plot(x, f(x))
% grid


figure
plot(x,y1)
hold on
plot(x,y2)

function y = phi(x,r_0,a)
y = r_0*tanh(a*x/r_0).*(x<=0) + (2-r_0)*tanh(a*x/(2-r_0)).*(x>0);
end

function y = sincy(x)
y = cos(x);
end