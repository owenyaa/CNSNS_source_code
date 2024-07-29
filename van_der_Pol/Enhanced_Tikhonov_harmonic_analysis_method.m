clear;clc;%close all;
tic;
global Tdata parameter_a N_dof N_harm N_w0 index_global N_min alpha temp_zeros
%% different initial values
% 不同系统系需要更改的参数
N_dof=1;N_harm=50;N_w0=2;
N_min=20; alpha=0.09;
% N_min=5; alpha=0.09;
% N_min=5; alpha=0.03;
Tdata=0:0.1:500;
%% 计算基频的组合系数
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
%index_global=[index_global,index_global,index_global];

%% 第一行存储频率，后面存储谐波系数,每个自由度两列
% parameter_a=[C_11,S_11,C_21,S_21;
%             C_12,S_12,C_22,S_22;
%             C_13,S_13,C_23,S_23;...];
parameter_a=zeros(N_harm,2*N_dof);
parameter_a(1:2,:)=[0.1,0.1;0.1,0.1];
%%
ini_parameter_a=parameter_a;
iteration=length(Tdata);
% to indicate where a is in the parametric space
parameter_a_judge=@(parameter_a)(abs(parameter_a(1,:))<1);
gammaT=1.414;rhob=0.5; % parameter for trust-region algorithm
parameter_a_record=parameter_a(1:N_harm,1:2*N_dof); % To record the values of parameters during iteration, and for each iteration, a is recorded in a single row of a_record
TR_record=[];  % recording the parameters during trust region
%% 本算例采用速度和加速度响应识别，且三个自由度噪声等级不同，线加速度1%，角加速度2%和5%
%% Response sensitivity iteration
Nmax=1000;   % maximum number for response sensitivity iteration
Ntr=20;      % maximum number for trust region iteration
%% response sensitivity Solution by ode45
for iii=1:Nmax
    %% 重新载入w0,此处为一个周期内均等选取1K个点
    % compute response and response sensitivity for each incremental
    Etol=1e-10;  % Relative error tolerance for convergence of the RS algorithm
    %%
    residual_iden=cal_residual(parameter_a);
    %% SSS为位移响应灵敏度矩阵，第一列和第二列为残差响应，频率灵敏度从第三列开始
    % 计算的灵敏度矩阵中包含S_11，但是这里要去掉
    SSS=reshape(residual_iden(:,N_dof+1:2*N_dof),N_dof*length(Tdata),1);%w01
    for i=1:2*N_harm*N_dof-1
        SSS=[SSS,reshape(residual_iden(:,N_dof*(i+1)+1:N_dof*(i+2)),N_dof*length(Tdata),1)];
    end
    
    dR=-reshape(residual_iden(:,1:N_dof),N_dof*length(Tdata),1);
    [U,s,V]=csvd(SSS);
    lambda_inverse=l_curve(U,s,dR);
    atemp=parameter_a;
    % trust-region algorithm
    for trust=1:Ntr
        %% 计算的da排序，为w0,C_11,S_11(缺,不参与迭代),C_12,S_12,按照参数矩阵的形式排列da
        real_da=tikhonov(U,s,V,dR,lambda_inverse);
        da=reshape(real_da,2,N_dof*N_harm);da=da';
        
        %%
        if ~parameter_a_judge(atemp+da)      % if updated a is not out of the parametric space, then, lambda should be increased until updated a is in ...
            lambda_inverse=lambda_inverse*gammaT; %  update of lambda
            continue;
        end
        %% 计算的da排序，为w0,C_11,S_11(缺,不参与迭代),C_12,S_12,按照参数矩阵的形式排列da
        real_da=tikhonov(U,s,V,dR,lambda_inverse);
        da=reshape(real_da,2,N_dof*N_harm);da=da';
        parameter_a=atemp+da;

        N_real_harm=ceil(N_min+(N_harm-N_min)*exp(-alpha*(iii-1)));
        Harm_parameter_a=parameter_a;
        %% 计算每一阶频率的振幅
        for j=1:N_dof
            for i=1:N_harm
                vector_amplitude(i,j)=sqrt(Harm_parameter_a(i,2*j-1)^2+Harm_parameter_a(i,2*j)^2);
            end
        end
        % 对振幅进行降序排列，并记录降序之前的位置
        % descend降序，ascend升序，position原来为位置，
        [descend_vector_amplitude,position]=sort(vector_amplitude,1,'descend');
        temp_zeros=zeros(size(Harm_parameter_a));
        for j=1:N_dof
            for i=1:N_real_harm
                temp_zeros(position(i,j),2*j-1)=1;
                temp_zeros(position(i,j),2*j)=1;
            end
        end
        parameter_a=temp_zeros.*parameter_a;
        
        da=temp_zeros.*da;
        ini_da=da(1:N_harm,1:2);
        for num_dof=1:N_dof-1
            ini_da=[ini_da;da(1:N_harm,2*num_dof+1:2*(num_dof+1))];%da(N_harm+1:2*N_harm,1:2);
        end
        ini_da=ini_da';ini_da=reshape(ini_da,2*N_dof*N_harm,1);

        %% 用新的parameter_a=atemp+da计算响应
        residual_da=cal_residual(parameter_a);
        dRtemp=-reshape(residual_da(:,1:N_dof),N_dof*length(Tdata),1);
        LdR=SSS*ini_da-dR;
        rhos=(dR'*dR-dRtemp'*dRtemp)/(dR'*dR-LdR'*LdR);  % agreement indicator
        if rhos>=rhob
            break;
        end
        lambda_inverse=lambda_inverse*gammaT;
    end
    tolt=norm(da)/norm(parameter_a);
    parameter_a_record=[parameter_a_record,parameter_a];
    TR_record=[TR_record;lambda_inverse];
    parameter_a
    if tolt<=Etol
        break;
    end
    every_a(iii).parameter_a=parameter_a;
    iii
end
toc;
residua=cal_residual(parameter_a);

figure;
plot(Tdata,residua(:,1),'r-','LineWidth',1);
set(gca,'FontName','Times New Roman','FontSize',15,'LineWidth',1.5);

% load 'N_9_50_20.mat';
Tdata=0:0.01:200;
x=zeros(N_dof,length(Tdata));dx=zeros(N_dof,length(Tdata));ddx=zeros(N_dof,length(Tdata));
% 有X,Y两个自由度
w1=1;w2=0.765;
fundamental_w=[w1,w2];
vector_w=index_global(:,1:N_w0)*fundamental_w';
Harm_parameter_a=parameter_a;
for j=1:N_dof
    for i=1:N_harm   % i=1,3,5
        x(j,:)=x(j,:)+Harm_parameter_a(i,2*j-1)*cos(vector_w(i)*Tdata)+Harm_parameter_a(i,2*j)*sin(vector_w(i)*Tdata);
        dx(j,:)=dx(j,:)-vector_w(i)*Harm_parameter_a(i,2*j-1)*sin(vector_w(i)*Tdata)+vector_w(i)*Harm_parameter_a(i,2*j)*cos(vector_w(i)*Tdata);
        ddx(j,:)=ddx(j,:)-(vector_w(i))^2*Harm_parameter_a(i,2*j-1)*cos(vector_w(i)*Tdata)-(vector_w(i))^2*Harm_parameter_a(i,2*j)*sin(vector_w(i)*Tdata);
    end
end
figure;
plot(Tdata,x(1,:),'r-','LineWidth',1);
set(gca,'FontName','Times New Roman','FontSize',15,'LineWidth',1.5);
figure;
plot(x(1,:),dx(1,:),'r-','LineWidth',1);
set(gca,'FontName','Times New Roman','FontSize',15,'LineWidth',1.5);

ini_parameter_a
% 
% for j=1:N_dof
%     for i=1:N_harm
%         vector_amplitude(i,j)=sqrt(Harm_parameter_a(i,2*j-1)^2+Harm_parameter_a(i,2*j)^2);
%     end
% end
% figure;
% plot(abs(vector_w),vector_amplitude,'k.','LineWidth',1.5);
% set(gca,'FontName','Times New Roman','FontSize',15,'LineWidth',1.5);














