%% Author : GUANG_LIU  * owenyaa@gmail.com *
% Created Time : 2016-11-03 17:10
% Last Revised : GUANG_LIU ,2016-11-03
% Remark : ������Ľṹ�����ܺͱ�����һ��˵��%

clear;clc;close all;
%% 
global A ep f1 k_n w M
ep=0.1;lamda=0.4;k_n=5;f1=0.3;w=1;

M=[1,0;0,ep];
C=[ep*lamda,-ep*lamda;-ep*lamda,ep*lamda];K=[1,0;0,0];
A=[zeros(2),eye(2);-M\K,-M\C];
% odex=[x(1,1);x(2,1);dx(1,1);dx(2,1)];
odex=[0.3;0;0;0];
tt=0:0.01:1000;
options=odeset('RelTol',1e-10,'AbsTol',1e-10);
[t,num]=ode45('odehomo',tt,odex,options);
% [t,num]=ode45('odehomo',tt,odex);

figure;
plot(tt,num(:,1),'k-');
hold on
plot(tt,num(:,2),'r-');
set(gca,'FontName','Times New Roman','FontSize',15,'LineWidth',1.5);



figure;
plot(tt(1:50:end),num(1:50:end,1),'k.');
hold on
plot(tt(1:50:end),num(1:50:end,2),'r.');
set(gca,'FontName','Times New Roman','FontSize',15,'LineWidth',1.5);
figure;
plot(num(70000:end,1),num(70000:end,3),'k-');
hold on
plot(num(70000:end,2),num(70000:end,4),'r-');
set(gca,'FontName','Times New Roman','FontSize',15,'LineWidth',1.5);


% residual=-M\(C*num(:,3:4)'+K*num(:,1:2)'+[-ep*f1*cos(w*tt)+ep*k_n*(num(:,1)'-num(:,2)').^3;ep*k_n*(num(:,2)'-num(:,1)').^3]);
% figure;
% plot(tt,residual(1,:)-ddx(1,:),'k-');
% hold on
% plot(tt,residual(2,:)-ddx(2,:),'r-');
% 
% figure;
% plot(tt,ddx(1,:),'r-');
% hold on
% plot(tt,ddx(2,:),'k-');

% dt=tt(2)-tt(1);
% data=num(1:end,1);
% N=length(data);
% N_fft=2^16;
% Y=fft(data,N_fft);
% Pyy=Y.*conj(Y)/N_fft;
% f=1/dt*(0:N_fft/2)/N_fft;
% figure��
% plot(f,Pyy(1:(N_fft/2+1))/10000,'r-','LineWidth',1);
% set(gca,'FontName','Times New Roman','FontSize',15,'LineWidth',1.5);


fs=100;%����Ƶ��
% ����Ƶ����ʱ����֮��Ĺ�ϵ�� fs=1/dt
% ��������������ǣ�����Ƶ��Ҫ�����ź�Ƶ�ʵ������� 
N=2^16;  %��������2^17
% N�������㣬����FFT֮�󣬾Ϳ��Եõ�N�����FFT�����Ϊ�˷������FFT���㣬ͨ��Nȡ2�������η���
% Ҫ��ȷ��xHz������Ҫ��������Ϊ1/x����źţ�����FFT��
% Ҫ���Ƶ�ʷֱ��ʣ�����Ҫ���Ӳ�������
n=0:N-1;
t=n/fs;  % dt=1/fs ��ʾʱ����   fs=1/dt
y=fft(num(700000:end,1),N);  % ����fft�任
% �������Ƶ��ΪFs���ź�Ƶ��F����������ΪN����ôFFT֮��������һ��ΪN��ĸ�����
% ÿһ����Ͷ�Ӧ��һ��Ƶ�ʵ㡣������ģֵ�����Ǹ�Ƶ��ֵ�µķ������ԡ�
% y % ���y����fft֮��Ľ����
m=abs(y(1:N/2))*2/N; % ���źŵ���ʵ��ֵ
f=2*pi*n*fs/N;  % ע�� ��2*pi �Ͳ���2*pi������%m=log10(m);
figure;
% plot(f(1:N/2)./(2*pi),m(1:N/2),'r-','LineWidth',1);
plot(f(1:N/2),m(1:N/2),'r-','LineWidth',1);
set(gca,'FontName','Times New Roman','FontSize',15,'LineWidth',1.5);
h1=legend('$$Iteration steps$$');
set(h1,'Interpreter','latex','FontSize',15);
set(gca,'FontName','Times New Roman','FontSize',15,'LineWidth',1.5);

y=fft(num(700000:end,2),N);  % ����fft�任
% �������Ƶ��ΪFs���ź�Ƶ��F����������ΪN����ôFFT֮��������һ��ΪN��ĸ�����
% ÿһ����Ͷ�Ӧ��һ��Ƶ�ʵ㡣������ģֵ�����Ǹ�Ƶ��ֵ�µķ������ԡ�
% y % ���y����fft֮��Ľ����
m=abs(y(1:N/2))*2/N; % ���źŵ���ʵ��ֵ
f=2*pi*n*fs/N;  % ע�� ��2*pi �Ͳ���2*pi������%m=log10(m);
figure;
% plot(f(1:N/2)./(2*pi),m(1:N/2),'r-','LineWidth',1);
plot(f(1:N/2),m(1:N/2),'k-','LineWidth',1);
set(gca,'FontName','Times New Roman','FontSize',15,'LineWidth',1.5);
h1=legend('$$Iteration steps$$');
set(h1,'Interpreter','latex','FontSize',15);
set(gca,'FontName','Times New Roman','FontSize',15,'LineWidth',1.5);

% hold on;
% plot(tt(1:30:end),num(1:30:end,1),'k.','MarkerSize',15);
% set(gca,'FontName','Times New Roman','FontSize',15,'LineWidth',1.5);
% figure;
% plot(num(1:30:end,1),num(1:30:end,2),'k.','MarkerSize',15);
% set(gca,'FontName','Times New Roman','FontSize',15,'LineWidth',1.5);


% aa=ddx+e*(1-x.^2).*dx+w0_2*x+a*x.^3-f1*cos(w1*Tdata)-f2*cos(w2*Tdata);
% figure;
% plot(num(1:end,1),num(1:end,2),'k-','LineWidth',1.5);
% set(gca,'FontName','Times New Roman','FontSize',15,'LineWidth',1.5);
