function residual=cal_residual(parameter_a)
global index_global N_dof Tdata N_harm N_w0 w
Harm_parameter_a=parameter_a(2:end,:);
%% 自由度数目，谐波数目，参数个数2*N_harm*N_dof
%系统参数
ep=0.1;lamda=0.4;k_n=5;f1=0.3;

M=[1,0;0,ep];C=[ep*lamda,-ep*lamda;-ep*lamda,ep*lamda];K=[1,0;0,0];

%% 计算两个基频率的组合，注意，频率组合要去掉负数频率
fundamental_w=parameter_a(1,1:N_w0);fundamental_w=[fundamental_w,fundamental_w];
for i=1:N_dof
    vector_w(:,i)=index_global(:,2*i-1:2*i)*fundamental_w(1,2*i-1:2*i)';
end
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
residual(1:N_dof,:)=M*ddx+C*dx+K*x+[-ep*f1*cos(w*Tdata)+ep*k_n*(x(1,:)-x(2,:)).^3;ep*k_n*(x(2,:)-x(1,:)).^3];
%% 计算频率的灵敏度  只需要计算基频的就可以
for j=1:N_w0-1
    x_w0=zeros(N_dof,length(Tdata));dx_w0=zeros(N_dof,length(Tdata));ddx_w0=zeros(N_dof,length(Tdata));
    for ij=1:N_dof
        for i=1:N_harm   % i=1,3,5
            x_w0(ij,:)=x_w0(ij,:)-index_global(i,j)*Harm_parameter_a(i,2*ij-1)*Tdata.*sin(vector_w(i,ij)*Tdata)+index_global(i,j)*Harm_parameter_a(i,2*ij)*Tdata.*cos(vector_w(i,ij)*Tdata);
            dx_w0(ij,:)=dx_w0(ij,:)-(index_global(i,j)*Harm_parameter_a(i,2*ij-1)*sin(vector_w(i,ij)*Tdata)+vector_w(i,ij)*index_global(i,j)*Tdata*Harm_parameter_a(i,2*ij-1).*cos(vector_w(i,ij)*Tdata))+...
                (index_global(i,j)*Harm_parameter_a(i,2*ij)*cos(vector_w(i,ij)*Tdata)-vector_w(i,ij)*index_global(i,j)*Tdata*Harm_parameter_a(i,2*ij).*sin(vector_w(i,ij)*Tdata));
            ddx_w0(ij,:)=ddx_w0(ij,:)-(2*vector_w(i,ij)*index_global(i,j)*Harm_parameter_a(i,2*ij-1)*cos(vector_w(i,ij)*Tdata)-vector_w(i,ij)^2*index_global(i,j)*Tdata*Harm_parameter_a(i,2*ij-1).*sin(vector_w(i,ij)*Tdata))-...
                (2*vector_w(i,ij)*index_global(i,j)*Harm_parameter_a(i,2*ij)*sin(vector_w(i,ij)*Tdata)+vector_w(i,ij)^2*index_global(i,j)*Tdata*Harm_parameter_a(i,2*ij).*cos(vector_w(i,ij)*Tdata));
        end
    end
    residual(N_dof*j+1:(1+j)*N_dof,:)=M*ddx_w0+C*dx_w0+K*x_w0+3*ep*k_n*[(x(1,:)-x(2,:)).^2.*(x_w0(1,:)-x_w0(2,:));(x(2,:)-x(1,:)).^2.*(x_w0(2,:)-x_w0(1,:))];
end
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
            x_a(k,:)=x_a(k,:)+sensitivity_parameter_a(j,2*k-1)*cos(vector_w(j,k)*Tdata)+sensitivity_parameter_a(j,2*k)*sin(vector_w(j,k)*Tdata);
            dx_a(k,:)=dx_a(k,:)-vector_w(j,k)*sensitivity_parameter_a(j,2*k-1)*sin(vector_w(j,k)*Tdata)+vector_w(j,k)*sensitivity_parameter_a(j,2*k)*cos(vector_w(j,k)*Tdata);
            ddx_a(k,:)=ddx_a(k,:)-(vector_w(j,k))^2*sensitivity_parameter_a(j,2*k-1)*cos(vector_w(j,k)*Tdata)-(vector_w(j,k))^2*sensitivity_parameter_a(j,2*k)*sin(vector_w(j,k)*Tdata);
        end
    end
    residual(N_dof*(i+N_w0)-1:N_dof*(i+N_w0),:)=M*ddx_a+C*dx_a+K*x_a+3*ep*k_n*[(x(1,:)-x(2,:)).^2.*(x_a(1,:)-x_a(2,:));(x(2,:)-x(1,:)).^2.*(x_a(2,:)-x_a(1,:))];
end

residual=residual';

