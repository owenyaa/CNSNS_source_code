clear;clc;close all;
e=1.08;w0_2=1;a=1;f1=0.2;f2=1;w1=1;w2=0.765;
load 'N_9_50_20_a_0.03.mat';
%% 计算两个基频率的组合
fundamental_w=[w1,w2];%fundamental_w=[fundamental_w,fundamental_w,fundamental_w];
for i=1:N_dof
    vector_w(:,i)=index_global(:,2*i-1:2*i)*fundamental_w(1,2*i-1:2*i)';
end
Tdata=0:0.01:200;

temp_b10=-1:0.01:1;
temp_b20=-1:0.01:1;
for ii=1:length(temp_b10)
    for jj=1:length(temp_b20)
        Harm_parameter_a=parameter_a;
        Harm_parameter_a(1,1)=temp_b10(ii);
        Harm_parameter_a(2,1)=temp_b20(jj);
        %% 计算方程残差
        x=zeros(N_dof,length(Tdata));dx=zeros(N_dof,length(Tdata));ddx=zeros(N_dof,length(Tdata));
        % 有三个自由度
        for j=1:N_dof
            for i=1:N_harm   % i=1,3,5
                x(j,:)=x(j,:)+Harm_parameter_a(i,2*j-1)*cos(vector_w(i,j)*Tdata)+Harm_parameter_a(i,2*j)*sin(vector_w(i,j)*Tdata);
                dx(j,:)=dx(j,:)-vector_w(i,j)*Harm_parameter_a(i,2*j-1)*sin(vector_w(i,j)*Tdata)+vector_w(i,j)*Harm_parameter_a(i,2*j)*cos(vector_w(i,j)*Tdata);
                ddx(j,:)=ddx(j,:)-(vector_w(i,j))^2*Harm_parameter_a(i,2*j-1)*cos(vector_w(i,j)*Tdata)-(vector_w(i,j))^2*Harm_parameter_a(i,2*j)*sin(vector_w(i,j)*Tdata);
            end
        end
        residual(1:N_dof,:)=ddx+e*(1-x.^2).*dx+w0_2*x+a*x.^3-f1*cos(w1*Tdata)-f2*cos(w2*Tdata);
        R10(ii,jj)=norm(residual);
    end
    ii
end
figure;
hold on
surf(temp_b10,temp_b20,-100+0*R10,R10);
h=surfc(temp_b10,temp_b20,R10);%,'r-','LineWidth',1);
h1=legend('$$a_0$$');
set(h1,'Interpreter','latex','FontSize',15);
set(gca,'FontName','Times New Roman','FontSize',15,'LineWidth',1.5);


[C1,h] = contour3(temp_b10,temp_b20,R10);
set(h,'ShowText','on','TextStep',get(h,'LevelStep')*0.2)
colormap cool
h1=legend('$$b_1$$','$$b_2$$');
set(h1,'Interpreter','latex','FontSize',15);


Harm_parameter_a=parameter_a;
% Harm_parameter_a(1,1)=temp_b10(ii);
% Harm_parameter_a(2,1)=temp_b20(jj);
%% 计算方程残差
x=zeros(N_dof,length(Tdata));dx=zeros(N_dof,length(Tdata));ddx=zeros(N_dof,length(Tdata));
% 有三个自由度
for j=1:N_dof
    for i=1:N_harm   % i=1,3,5
        x(j,:)=x(j,:)+Harm_parameter_a(i,2*j-1)*cos(vector_w(i,j)*Tdata)+Harm_parameter_a(i,2*j)*sin(vector_w(i,j)*Tdata);
        dx(j,:)=dx(j,:)-vector_w(i,j)*Harm_parameter_a(i,2*j-1)*sin(vector_w(i,j)*Tdata)+vector_w(i,j)*Harm_parameter_a(i,2*j)*cos(vector_w(i,j)*Tdata);
        ddx(j,:)=ddx(j,:)-(vector_w(i,j))^2*Harm_parameter_a(i,2*j-1)*cos(vector_w(i,j)*Tdata)-(vector_w(i,j))^2*Harm_parameter_a(i,2*j)*sin(vector_w(i,j)*Tdata);
    end
end
residual(1:N_dof,:)=ddx+e*(1-x.^2).*dx+w0_2*x+a*x.^3-f1*cos(w1*Tdata)-f2*cos(w2*Tdata);
R110=norm(residual);
hold on
plot3(Harm_parameter_a(2,1),Harm_parameter_a(1,1),R110,'k.','MarkerSize',15);



