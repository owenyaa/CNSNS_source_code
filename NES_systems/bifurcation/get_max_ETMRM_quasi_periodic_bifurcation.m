clear;%close all;
load 'new_W_0945_1055_quasi.mat';

N_dof=2;N_harm_quasi=72;N_w0=2;%基频个数
N_min=35; alpha=0.02;

index_global=[1 0;0 1];
for N=3:2:11
% for N=2:1:5
    temp=[N 0];
    for n=N-1:-1:1
        m=N-n;
        temp=[temp;n m;n -m];
    end
    temp=[temp;0 N];
    index_global=[index_global;temp];
end
index_global=[index_global,index_global];
vector_w_w=0.945:0.005:1.055;

load 'W_08_12_period.mat';
% global parameter_a w N_harm N_dof
N_harm=10;

for ii=1:length(vector_w_w)
    w=every_a1(ii).w;
    parameter_a=every_a1(ii).parameter_a;
    Tdata=0:0.01:500;
    x=zeros(N_dof,length(Tdata));dx=zeros(N_dof,length(Tdata));ddx=zeros(N_dof,length(Tdata));
    % 有X,Y两个自由度
    fundamental_w=parameter_a(1,1:N_w0);
    vector_w=index_global(:,1:N_w0)*fundamental_w';
    Harm_parameter_a=parameter_a(2:end,:);
    for j=1:N_dof
        for i=1:N_harm_quasi   % i=1,3,5
            %if vector_w(i)>0
            x(j,:)=x(j,:)+Harm_parameter_a(i,2*j-1)*cos(vector_w(i)*Tdata)+Harm_parameter_a(i,2*j)*sin(vector_w(i)*Tdata);
            dx(j,:)=dx(j,:)-vector_w(i)*Harm_parameter_a(i,2*j-1)*sin(vector_w(i)*Tdata)+vector_w(i)*Harm_parameter_a(i,2*j)*cos(vector_w(i)*Tdata);
            ddx(j,:)=ddx(j,:)-(vector_w(i))^2*Harm_parameter_a(i,2*j-1)*cos(vector_w(i)*Tdata)-(vector_w(i))^2*Harm_parameter_a(i,2*j)*sin(vector_w(i)*Tdata);
            %end
        end
    end
    x=x';
    dx=dx';
    
    
    xmax_1=x(:,1);
    xmax_2=x(:,2);
    
    xmax1(ii)=max(getmax(xmax_1));
    xmax2(ii)=max(getmax(xmax_2));
    %     %% plot x_1
    %     xmax(ii).xmax=getmax(xmax_1);
    %     hold on;
    %     ww=w*ones(1,length(xmax(ii).xmax));
    %     plot(ww,xmax(ii).xmax,'k.','MarkerSize',8);
    %     hold on;
    %% plot x_2
    %     xmax(ii).xmax=getmax(xmax_2);
    %     hold on;
    %     ww=w*ones(1,length(xmax(ii).xmax));
    %     plot(ww,xmax(ii).xmax,'r.','MarkerSize',8);
    %     hold on;
end
set(gca,'FontName','Times New Roman','FontSize',15,'LineWidth',1.5);

for jj=1:81
    %peridoc
    w0=every_a(jj).w;parameter_a_peridoc=every_a(jj).parameter_a;
    w_peridoc(jj)=every_a(jj).w;
    Harm_parameter_a_peridoc=parameter_a_peridoc;
    %% 计算方程残差
    x_peridoc=zeros(N_dof,length(Tdata));dx_peridoc=zeros(N_dof,length(Tdata));ddx_peridoc=zeros(N_dof,length(Tdata));
    % 有X,Y两个自由度
    % for i=1:N_harm   % i=1,3,5
    for j=1:N_dof
        for i=1:N_harm   % i=1,3,5
            x_peridoc(j,:)=x_peridoc(j,:)+Harm_parameter_a_peridoc(i,2*j-1)*cos((2*i-1)*w0*Tdata)+Harm_parameter_a_peridoc(i,2*j)*sin((2*i-1)*w0*Tdata);
            dx_peridoc(j,:)=dx_peridoc(j,:)-w0*(2*i-1)*Harm_parameter_a_peridoc(i,2*j-1)*sin((2*i-1)*w0*Tdata)+w0*(2*i-1)*Harm_parameter_a_peridoc(i,2*j)*cos((2*i-1)*w0*Tdata);
            ddx_peridoc(j,:)=ddx_peridoc(j,:)-(w0*(2*i-1))^2*Harm_parameter_a_peridoc(i,2*j-1)*cos((2*i-1)*w0*Tdata)-(w0*(2*i-1))^2*Harm_parameter_a_peridoc(i,2*j)*sin((2*i-1)*w0*Tdata);
        end
    end
    x_peridoc=x_peridoc';
    dx_peridoc=dx_peridoc';
    
    
    x_peridocmax_1=x_peridoc(:,1);
    x_peridocmax_2=x_peridoc(:,2);
    
    x_peridocmax1(jj)=max(getmax(x_peridocmax_1));
    x_peridocmax2(jj)=max(getmax(x_peridocmax_2));
end


figure;
plot(w_peridoc(1,1:30),x_peridocmax1(1,1:30),'k-','LineWidth',1.5);
hold on
plot(w_peridoc(1,30:52),x_peridocmax1(1,30:52),'k--','LineWidth',1.5);
hold on
plot(w_peridoc(1,52:end),x_peridocmax1(1,52:end),'k-','LineWidth',1.5);
hold on
plot(vector_w_w,xmax1,'k-','LineWidth',1);
h1=legend('$$stable$$','$$unstable$$');
set(h1,'Interpreter','latex','FontSize',15);
set(gca,'FontName','Times New Roman','FontSize',15,'LineWidth',1.5);
figure;
plot(w_peridoc(1,1:30),x_peridocmax2(1,1:30),'r-','LineWidth',1.5);
hold on
plot(w_peridoc(1,30:52),x_peridocmax2(1,30:52),'r--','LineWidth',1.5);
hold on
plot(w_peridoc(1,52:end),x_peridocmax2(1,52:end),'r-','LineWidth',1.5);
hold on;
plot(vector_w_w,xmax2,'r-','LineWidth',1);
h1=legend('$$LC stable$$','$$LC unstable$$','$$QP stable$$');
set(h1,'Interpreter','latex','FontSize',15);
set(gca,'FontName','Times New Roman','FontSize',15,'LineWidth',1.5);
 
figure;
plot(Tdata,x(:,1),'k-','LineWidth',1);
hold on;
plot(Tdata,x(:,2),'b-','LineWidth',1);

figure;
plot(x(:,1),dx(:,1),'r-','LineWidth',1);
hold on;
plot(x(:,2),dx(:,2),'k-','LineWidth',1);
h1=legend('$$Extremes of x_1$$','$$Extremes of x_2$$');
set(h1,'Interpreter','latex','FontSize',15);
set(gca,'FontName','Times New Roman','FontSize',15,'LineWidth',1.5);








