clear;
%%        coefficient保存有所有谐波系数，对应风速和频率。coefficient.Q为风速，%
%          coefficient.A0_1为第一个自由度谐波系数，前5个位余弦系数，后5个位正弦系数。%
%          coefficient.A0_2和coefficient.A0_3同上。coefficient.frequency为对应的频率。%

% load coefficient_2_test.mat;
clear;close all;%clc;
load W_08_12_period.mat;
global A K M parameter_a w N_harm ep k_n f1 N_dof
N_harm=10;N_dof=2;
multipliers_mat=[];
ep=0.1;lamda=0.4;k_n=5;f1=0.3;%w=0.8;

M=[1,0;0,ep];C=[ep*lamda,-ep*lamda;-ep*lamda,ep*lamda];K=[1,0;0,0];
A=[zeros(2),eye(2);-M\K,-M\C];
for j=1:81
    %临界值 27超出(w=0.93)0.95 54进入(w=1.065)
    w(j)=every_a(j).w;parameter_a=every_a(j).parameter_a;
    %找出了系统的频率，并计算周期T,frequency为系统圆频率，不用再用2*pi除
    T=2*pi/w(j);
    A_11(j)=sqrt(parameter_a(1,1)^2+parameter_a(1,2)^2);
    A_21(j)=sqrt(parameter_a(1,3)^2+parameter_a(1,4)^2);
end
figure;
plot(w,A_11,'r-','LineWidth',1);
hold on;
plot(w,A_21,'k-','LineWidth',1);
h1=legend('$$x_1$$','$$x_2$$');
set(h1,'Interpreter','latex','FontSize',15);
set(gca,'FontName','Times New Roman','FontSize',15,'LineWidth',1.5);



t=0:0.01:2*pi;r=1+t.*0;
hold on;
polar(t,r);axis equal;
% for j=1:22
%     QQ(j)=coefficient(j).Q;
% end
% figure;
% plot(QQ,abs(multipliers_mat))

