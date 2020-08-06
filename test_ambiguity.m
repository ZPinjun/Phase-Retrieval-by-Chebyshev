clc;clear;close all
M = 50;
N = 100;

x = zeros(M,1);
x([1,40,69]) = [1.2;1.6;2.8];
y = abs(fft(x,N)).^2;
figure(1)
plot(1:N,y);hold on
x2 = zeros(M,1);
x2([1,30,69]) = [1.2;1.6;2.8];
y = abs(fft(x2,N)).^2;
plot(1:N,y);hold on