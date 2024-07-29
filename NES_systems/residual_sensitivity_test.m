clc;clear;close all;
global tf h N_dof N_harm N_w0 Tdata index_global w
%% parameters
w1=0.8533;
w=1;
N_dof=2;N_harm=32;N_w0=2;%基频个数
% 其他参数
w0=zeros(1,2*N_dof);
w0(1,1:N_w0)=[w1,w];
% 其他参数
tf=2*pi/min(w0(1,1:2));h=2*pi/(min(w0(1,1:2))*1000);Tdata=(0:h:tf);
%% 计算基频的组合系数
index_global=[1 0;0 1];
for N=3:2:7
    temp=[N 0];
    for n=N-1:-1:1
        m=N-n;
        temp=[temp;n m;n -m];
    end
    temp=[temp;0 N];
    index_global=[index_global;temp];
end
index_global=[index_global,index_global];
%% 第一行存储频率，后面存储谐波系数,每个自由度两列
% parameter_a=[w0,  0,    0,   0;
%             C_11,S_11,C_21,S_21;
%             C_12,S_12,C_22,S_22;
%             C_13,S_13,C_23,S_23;...];
% parameter_a(2,:)=[0.4,0.1,0.1,0.1,0.2,0.1,];% 一阶谐波初值

C0_1=[0.342172594555505,-0.000133005831228,0.000115963368886]';S0_1=[0.1,-0.000712982230221,-0.000108784126763]';
C0_2=[0.122467788388526,0.001563371451967,-0.000646804589735]';S0_2=[0.072813577124383,0.011674778196277,0.000571448020811]';
parameter_a=zeros(N_harm+1,2*N_dof);
parameter_a(1,:)=w0;
parameter_a(2:4,:)=[C0_1,S0_1,C0_2,S0_2];   % 6行2列

%% 计算残差
residual=cal_residual(parameter_a);
%% 绘制残差曲线
% figure;
% plot(Tdata,residual(:,1),'r-','LineWidth',1);
% hold on;
% plot(Tdata,residual(:,2),'k-','LineWidth',1);
% hold on;
% plot(Tdata,residual(:,3),'b-','LineWidth',1);
% h1=legend('$$h$$','$$\alpha$$','$$\beta$$');
% set(h1,'Interpreter','latex','FontSize',15);
% set(gca,'FontName','Times New Roman','FontSize',15,'LineWidth',1.5);

%% 此处为灵敏度验证程序，本质为用差分代替灵敏度，可以验证灵敏度是否算错
% 计算正问题之后，将某个参数减去一个小量，用新参数再算一次
% 两次所得结果做差再除以小量即为灵敏度，验证差分的灵敏度和直接计算的灵敏度曲线是否重合
% 但是要注意，差分结果肯定是对的(即下文差量结果)，可能出错的是原程序的内容，x_cal中内容
% 另外关于对应，差分结果(parameter_a1)中的位移(x1(1,:))，速度，加速度对应到原程序(parameter_a)中的灵敏度内容x_cal(2,:)
% 
for i=1:2*N_harm*N_dof
% for i=1:2
    ddt=0.000001;
    %初始化参数
    parameter_a=zeros(N_harm+1,2*N_dof);
    parameter_a(1,:)=w0;
    parameter_a(2:4,:)=[C0_1,S0_1,C0_2,S0_2];   % 6行2列
    
    %% 参数系数矩阵，1代表求此处参数的差分
    sensitivity_parameter_a1=zeros(2*N_harm*N_dof,1);
    sensitivity_parameter_a1(i,1)=1;
    sensitivity_parameter_a1=reshape(sensitivity_parameter_a1,2,N_harm*N_dof);sensitivity_parameter_a1=sensitivity_parameter_a1';
    sensitivity_parameter_a=sensitivity_parameter_a1(1:N_harm,1:2);
    for num_dof=1:N_dof-1
        sensitivity_parameter_a=[sensitivity_parameter_a,sensitivity_parameter_a1(num_dof*N_harm+1:(num_dof+1)*N_harm,1:2)];%da(N_harm+1:2*N_harm,1:2);
    end
    temp_real_w0=zeros(1,2*N_dof);
    sensitivity_parameter_a=[temp_real_w0;sensitivity_parameter_a];
    parameter_a=parameter_a+ddt*sensitivity_parameter_a;% 增量加上去了
    
    residual1=cal_residual(parameter_a);
    aaaa=(residual1-residual)/ddt;
    figure; 
    plot(Tdata,residual(:,N_dof*(i+N_w0)-1),'r-')
    hold on
    plot(Tdata,residual(:,N_dof*(i+N_w0)),'r-')
    hold on
    
    plot(Tdata,aaaa(:,1),'k-')
    hold on
    plot(Tdata,aaaa(:,2),'k-')
end


% for i=1:1
%     % i=8;
%     ddt=0.0001;
%     w0=zeros(1,2*N_dof);
%     w0(1,1:N_w0)=[w1,w];
%     w0(1,i)=w0(1,i)+ddt;
%     parameter_a(1,:)=w0;
%     parameter_a(2:4,:)=[C0_1,S0_1,C0_2,S0_2];   % 6行2列
%     
%     residual1=cal_residual(parameter_a);
%     aaaa=(residual1-residual)/ddt;
%     figure; 
%     plot(Tdata,residual(:,i*N_dof+1),'r-')
%     hold on
%     plot(Tdata,residual(:,i*N_dof+2),'r-')
%     hold on
%     
%     plot(Tdata,aaaa(:,1),'k-')
%     hold on
%     plot(Tdata,aaaa(:,2),'k-')
% end
















