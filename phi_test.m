clear all
close all


x = normrnd(0,1,50,1);
y = rand(50,1);

z = corr(x,y);