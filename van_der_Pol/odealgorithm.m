%% Author : GUANG_LIU  * owenyaa@gmail.com *
% Created Time : 2016-11-03 17:10
% Last Revised : GUANG_LIU ,2016-11-03
% Remark : 本程序的结构、功能和变量做一下说明%

% clear;clc;close all;
%% 
global A e a f1 f2 w1 w2 
e=1.08;w0_2=1;a=1;f1=0.2;f2=1;w1=1;w2=0.765;
A=[0,1;-w0_2,-e];

% odex=[x(1);dx(1)];
odex=[1.5;0];
tt=0:0.01:200;
options=odeset('RelTol',1e-6,'AbsTol',1e-6);
[t,num]=ode45('odehomo',tt,odex,options);
figure;
plot(tt,num(:,1),'k-','LineWidth',1.5);
h1=legend('$$x_1$$');
set(h1,'Interpreter','latex','FontSize',15);
set(gca,'FontName','Times New Roman','FontSize',15,'LineWidth',1.5);
figure;
plot(num(:,1),num(:,2),'k-','LineWidth',1.5);
h1=legend('$$x_1$$');
set(h1,'Interpreter','latex','FontSize',15);
set(gca,'FontName','Times New Roman','FontSize',15,'LineWidth',1.5);



% dt=tt(2)-tt(1);
% data=num(1:end,1);
% N=length(data);
% N_fft=2^16;
% Y=fft(data,N_fft);
% Pyy=Y.*conj(Y)/N_fft;
% f=1/dt*(0:N_fft/2)/N_fft;
% figure
% plot(f,Pyy(1:(N_fft/2+1))/10000,'r-','LineWidth',1);
% set(gca,'FontName','Times New Roman','FontSize',15,'LineWidth',1.5);


fs=100;%采样频率
% 采样频率与时间间隔之间的关系： fs=1/dt
% 采样定理告诉我们，采样频率要大于信号频率的两倍。 
N=2^16;  %采样点数2^17
% N个采样点，经过FFT之后，就可以得到N个点的FFT结果。为了方便进行FFT运算，通常N取2的整数次方。
% 要精确到xHz，则需要采样长度为1/x秒的信号，并做FFT。
% 要提高频率分辨率，就需要增加采样点数
n=0:N-1;
t=n/fs;  % dt=1/fs 表示时间间隔   fs=1/dt
y=fft(num(1:end,1),N);  % 进行fft变换
% 假设采样频率为Fs，信号频率F，采样点数为N。那么FFT之后结果就是一个为N点的复数。
% 每一个点就对应着一个频率点。这个点的模值，就是该频率值下的幅度特性。
% y % 输出y看看fft之后的结果。
m=abs(y(1:N/2))*2/N; % 求信号的真实幅值
f=2*pi*n*fs/N;  % 注意 乘2*pi 和不乘2*pi的区别%m=log10(m);
figure;
plot(f(1:N/2),m(1:N/2),'r-','LineWidth',1);
set(gca,'FontName','Times New Roman','FontSize',15,'LineWidth',1.5);
h1=legend('$$Iteration steps$$');
set(h1,'Interpreter','latex','FontSize',15);
set(gca,'FontName','Times New Roman','FontSize',15,'LineWidth',1.5);


% hold on;
plot(tt(1:30:end),num(1:30:end,1),'k.','MarkerSize',15);
set(gca,'FontName','Times New Roman','FontSize',15,'LineWidth',1.5);
figure;
plot(num(1:30:end,1),num(1:30:end,2),'k.','MarkerSize',15);
set(gca,'FontName','Times New Roman','FontSize',15,'LineWidth',1.5);


% aa=ddx+e*(1-x.^2).*dx+w0_2*x+a*x.^3-f1*cos(w1*Tdata)-f2*cos(w2*Tdata);
% figure;
% plot(num(1:end,1),num(1:end,2),'k-','LineWidth',1.5);
% set(gca,'FontName','Times New Roman','FontSize',15,'LineWidth',1.5);
