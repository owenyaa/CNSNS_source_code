clear;
%%        coefficient����������г��ϵ������Ӧ���ٺ�Ƶ�ʡ�coefficient.QΪ���٣�%
%          coefficient.A0_1Ϊ��һ�����ɶ�г��ϵ����ǰ5��λ����ϵ������5��λ����ϵ����%
%          coefficient.A0_2��coefficient.A0_3ͬ�ϡ�coefficient.frequencyΪ��Ӧ��Ƶ�ʡ�%

% load coefficient_2_test.mat;
clear;close all;%clc;
load period_solution_0.8_1.2.mat;
global A K M parameter_a w N_harm ep k_n f1 N_dof
N_harm=10;N_dof=2;
multipliers_mat=[];
ep=0.1;lamda=0.4;k_n=5;f1=0.3;%w=0.8;

M=[1,0;0,ep];C=[ep*lamda,-ep*lamda;-ep*lamda,ep*lamda];K=[1,0;0,0];
A=[zeros(2),eye(2);-M\K,-M\C];
for j=1:81
    %�ٽ�ֵ 27����(w=0.93)0.95 54����(w=1.065)
    w=every_a(j).w;parameter_a=every_a(j).parameter_a;
    %�ҳ���ϵͳ��Ƶ�ʣ�����������T,frequencyΪϵͳԲƵ�ʣ���������2*pi��
    T=2*pi/w;
    dx0=zeros(4,4);dx_t=zeros(4,4);
    for i=1:4
        dx0(i,i)=1;
    end
    tt=0:T/10000:T;
    for i=1:4
        [t,num]=ode45('ode_floquet',tt,dx0(:,i));
        temp_dx=num(end,:);
        dx_t(:,i)=temp_dx';
    end
    %temp_num(j,:)=num(length(num(:,1)),:);
    %D_tr_dx=eig(dx_t);
    floquet_multipliers=eig(dx_t)
    plot(real(floquet_multipliers),imag(floquet_multipliers),'k.')
    hold on
    
    multipliers_mat=[multipliers_mat floquet_multipliers];
    j
end
t=0:0.01:2*pi;r=1+t.*0;
hold on;
polar(t,r);axis equal;
% for j=1:22
%     QQ(j)=coefficient(j).Q;
% end
% figure;
% plot(QQ,abs(multipliers_mat))

