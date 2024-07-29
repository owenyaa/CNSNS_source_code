clear;clc;%close all;
%% different initial values
% 不同系统系需要更改的参数
global Tdata parameter_a %N_dof N_harm N_w0 index_global N_min alpha
%% different initial values
% 不同系统系需要更改的参数
% N_dof=1;N_harm=32;N_w0=20;N_min=10; alpha=0.03;
% Tdata=0:0.1:500;
% %% 计算基频的组合系数
index_global=[1 0;0 1];
for N=3:2:9
    temp=[N 0];
    for n=N-1:-1:1
        m=N-n;
        temp=[temp;n m;n -m];
    end
    temp=[temp;0 N];
    index_global=[index_global;temp];
end

% load N_3_8_6.mat;
% residua=cal_residual(parameter_a);
% 
% subplot(4,1,1);
% plot(Tdata,residua(:,1),'k-','LineWidth',1);
% set(gca,'FontName','Times New Roman','FontSize',15,'LineWidth',1.5);
figure;
load N_9_50_20_a_0.03.mat;
residua=cal_residual(parameter_a);

subplot(3,1,1);
plot(Tdata,residua(:,1),'k-','LineWidth',1);
set(gca,'FontName','Times New Roman','FontSize',15,'LineWidth',1.5);

load N_9_50_20_a_0.09.mat;
residua=cal_residual(parameter_a);

subplot(3,1,2);
plot(Tdata,residua(:,1),'k-','LineWidth',1);
set(gca,'FontName','Times New Roman','FontSize',15,'LineWidth',1.5);

load N_9_50_5_a_0.09.mat;
residua=cal_residual(parameter_a);

% load 'N_9_50_20.mat';
% N_dof=1;N_w0=2;
% Tdata=0:0.01:200;
% x=zeros(N_dof,length(Tdata));dx=zeros(N_dof,length(Tdata));ddx=zeros(N_dof,length(Tdata));
% % 有X,Y两个自由度
% w1=1;w2=0.765;
% fundamental_w=[w1,w2];
% vector_w=index_global(:,1:N_w0)*fundamental_w';
% Harm_parameter_a=parameter_a;
% for j=1:N_dof
%     for i=1:N_harm   % i=1,3,5
%         x(j,:)=x(j,:)+Harm_parameter_a(i,2*j-1)*cos(vector_w(i)*Tdata)+Harm_parameter_a(i,2*j)*sin(vector_w(i)*Tdata);
%         dx(j,:)=dx(j,:)-vector_w(i)*Harm_parameter_a(i,2*j-1)*sin(vector_w(i)*Tdata)+vector_w(i)*Harm_parameter_a(i,2*j)*cos(vector_w(i)*Tdata);
%         ddx(j,:)=ddx(j,:)-(vector_w(i))^2*Harm_parameter_a(i,2*j-1)*cos(vector_w(i)*Tdata)-(vector_w(i))^2*Harm_parameter_a(i,2*j)*sin(vector_w(i)*Tdata);
%     end
% end
% figure;
% plot(Tdata,x(1,:),'r-','LineWidth',1);
% set(gca,'FontName','Times New Roman','FontSize',15,'LineWidth',1.5);
% figure;
% plot(x(1,:),dx(1,:),'r-','LineWidth',1);
% set(gca,'FontName','Times New Roman','FontSize',15,'LineWidth',1.5);
% 
% ini_parameter_a
% 
% for j=1:N_dof
%     for i=1:N_harm
%         vector_amplitude(i,j)=sqrt(Harm_parameter_a(i,2*j-1)^2+Harm_parameter_a(i,2*j)^2);
%     end
% end
% figure;
% plot(abs(vector_w),vector_amplitude,'k.','LineWidth',1.5);
% set(gca,'FontName','Times New Roman','FontSize',15,'LineWidth',1.5);
% 





subplot(3,1,3);
plot(Tdata,residua(:,1),'k-','LineWidth',1);
set(gca,'FontName','Times New Roman','FontSize',15,'LineWidth',1.5);
