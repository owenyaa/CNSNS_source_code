function residual=cal_residual(parameter_a)
global N_dof Tdata N_harm w
Harm_parameter_a=parameter_a;
%% 自由度数目，谐波数目，参数个数2*N_harm*N_dof
%系统参数
ep=0.1;lamda=0.4;k_n=5;f1=0.3;

M=[1,0;0,ep];C=[ep*lamda,-ep*lamda;-ep*lamda,ep*lamda];K=[1,0;0,0];
%% 计算方程残差
x=zeros(N_dof,length(Tdata));dx=zeros(N_dof,length(Tdata));ddx=zeros(N_dof,length(Tdata));
% 有X,Y两个自由度
for j=1:N_dof
    for i=1:N_harm   % i=1,3,5
        x(j,:)=x(j,:)+Harm_parameter_a(i,2*j-1)*cos((2*i-1)*w*Tdata)+Harm_parameter_a(i,2*j)*sin((2*i-1)*w*Tdata);
        dx(j,:)=dx(j,:)-w*(2*i-1)*Harm_parameter_a(i,2*j-1)*sin((2*i-1)*w*Tdata)+w*(2*i-1)*Harm_parameter_a(i,2*j)*cos((2*i-1)*w*Tdata);
        ddx(j,:)=ddx(j,:)-(w*(2*i-1))^2*Harm_parameter_a(i,2*j-1)*cos((2*i-1)*w*Tdata)-(w*(2*i-1))^2*Harm_parameter_a(i,2*j)*sin((2*i-1)*w*Tdata);
    end
end
residual(1:N_dof,:)=M*ddx+C*dx+K*x+[-ep*f1*cos(w*Tdata)+ep*k_n*(x(1,:)-x(2,:)).^3;ep*k_n*(x(2,:)-x(1,:)).^3];
%% 计算频率的灵敏度  只需要计算基频的就可以
% x_w0=zeros(N_dof,length(Tdata));dx_w0=zeros(N_dof,length(Tdata));ddx_w0=zeros(N_dof,length(Tdata));
% for k=1:N_dof
%     for i=1:N_harm   % i=1,3,5
%         x_w0(k,:)=x_w0(k,:)-Harm_parameter_a(i,2*k-1)*(2*i-1)*Tdata.*sin((2*i-1)*w0*Tdata)+Harm_parameter_a(i,2*k)*(2*i-1)*Tdata.*cos((2*i-1)*w0*Tdata);
%         dx_w0(k,:)=dx_w0(k,:)-((2*i-1)*Harm_parameter_a(i,2*k-1)*sin((2*i-1)*w0*Tdata)+w0*Tdata*(2*i-1)^2*Harm_parameter_a(i,2*k-1).*cos((2*i-1)*w0*Tdata))+...
%             ((2*i-1)*Harm_parameter_a(i,2*k)*cos((2*i-1)*w0*Tdata)-w0*Tdata*(2*i-1)^2*Harm_parameter_a(i,2*k).*sin((2*i-1)*w0*Tdata));
%         ddx_w0(k,:)=ddx_w0(k,:)-(2*w0*(2*i-1)^2*Harm_parameter_a(i,2*k-1)*cos((2*i-1)*w0*Tdata)-w0^2*Tdata*(2*i-1)^3*Harm_parameter_a(i,2*k-1).*sin((2*i-1)*w0*Tdata))-...
%             (2*w0*(2*i-1)^2*Harm_parameter_a(i,2*k)*sin((2*i-1)*w0*Tdata)+w0^2*Tdata*(2*i-1)^3*Harm_parameter_a(i,2*k).*cos((2*i-1)*w0*Tdata));
%     end
% end
% residual(N_dof*j+1:(1+j)*N_dof,:)=M*ddx_w0+C*dx_w0+K*x_w0+3*ep*k_n*[(x(1,:)-x(2,:)).^2.*(x_w0(1,:)-x_w0(2,:));(x(2,:)-x(1,:)).^2.*(x_w0(2,:)-x_w0(1,:))];
%% 计算谐波系数的灵敏度
for i=1:2*N_harm*N_dof
    sensitivity_parameter_a1=zeros(2*N_harm*N_dof,1);
    x_a=zeros(N_dof,length(Tdata));dx_a=zeros(N_dof,length(Tdata));ddx_a=zeros(N_dof,length(Tdata));
    sensitivity_parameter_a1(i,1)=1;
    sensitivity_parameter_a1=reshape(sensitivity_parameter_a1,2,N_harm*N_dof);sensitivity_parameter_a1=sensitivity_parameter_a1';
    sensitivity_parameter_a=sensitivity_parameter_a1(1:N_harm,1:2);
    for num_dof=1:N_dof-1
        sensitivity_parameter_a=[sensitivity_parameter_a,sensitivity_parameter_a1(num_dof*N_harm+1:(num_dof+1)*N_harm,1:2)];%da(N_harm+1:2*N_harm,1:2);
    end
    for k=1:N_dof
        for j=1:N_harm
            x_a(k,:)=x_a(k,:)+sensitivity_parameter_a(j,2*k-1)*cos((2*j-1)*w*Tdata)+sensitivity_parameter_a(j,2*k)*sin((2*j-1)*w*Tdata);
            dx_a(k,:)=dx_a(k,:)-w*(2*j-1)*sensitivity_parameter_a(j,2*k-1)*sin((2*j-1)*w*Tdata)+w*(2*j-1)*sensitivity_parameter_a(j,2*k)*cos((2*j-1)*w*Tdata);
            ddx_a(k,:)=ddx_a(k,:)-(w*(2*j-1))^2*sensitivity_parameter_a(j,2*k-1)*cos((2*j-1)*w*Tdata)-(w*(2*j-1))^2*sensitivity_parameter_a(j,2*k)*sin((2*j-1)*w*Tdata);
        end
    end
    residual(N_dof*i+1:N_dof*(i+1),:)=M*ddx_a+C*dx_a+K*x_a+3*ep*k_n*[(x(1,:)-x(2,:)).^2.*(x_a(1,:)-x_a(2,:));(x(2,:)-x(1,:)).^2.*(x_a(2,:)-x_a(1,:))];
end

residual=residual';

