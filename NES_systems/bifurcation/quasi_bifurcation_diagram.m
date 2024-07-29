clear;close all;clc;
load new_W_0995_106_quasi.mat;
% global K M parameter_a w N_harm ep k_n f1 N_dof
% N_harm=10;N_dof=2;
% multipliers_mat=[];
% ep=0.1;lamda=0.4;k_n=5;f1=0.3;%w=0.8;
% 
% M=[1,0;0,ep];C=[ep*lamda,-ep*lamda;-ep*lamda,ep*lamda];K=[1,0;0,0];
w=0.95:0.005:1.06;
% w=every_a1.w;
for j=1:length(w)
    %临界值 27超出(w=0.93)0.95 54进入(w=1.065)
    parameter_a=every_a1(j).parameter_a;
    %找出了系统的频率，并计算周期T,frequency为系统圆频率，不用再用2*pi除
    %T=2*pi/w(j);
    %这个是外激励对应的振幅
    A_11(j)=sqrt(parameter_a(2,1)^2+parameter_a(2,2)^2);
    A_21(j)=sqrt(parameter_a(2,3)^2+parameter_a(2,4)^2);
end
figure;
plot(w,A_11,'k-','LineWidth',1);
hold on;
plot(w,A_21,'r-','LineWidth',1);
h1=legend('$$x_1$$','$$x_2$$');
set(h1,'Interpreter','latex','FontSize',15);
set(gca,'FontName','Times New Roman','FontSize',15,'LineWidth',1.5);
