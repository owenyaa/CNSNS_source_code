function residual=cal_residual(parameter_a)
global index_global N_dof Tdata N_harm
e=1.08;w0_2=1;a=1;f1=0.2;f2=1;w1=1;w2=0.765;
Harm_parameter_a=parameter_a;
%% 自由度数目，谐波数目，参数个数2*N_harm*N_dof

%% 计算两个基频率的组合
fundamental_w=[w1,w2];%fundamental_w=[fundamental_w,fundamental_w,fundamental_w];
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
residual(1:N_dof,:)=ddx+e*(1-x.^2).*dx+w0_2*x+a*x.^3-f1*cos(w1*Tdata)-f2*cos(w2*Tdata);
%% 计算频率的灵敏度  只需要计算基频的就可以
% for j=1:N_w0
%     x_w0=zeros(N_dof,length(Tdata));dx_w0=zeros(N_dof,length(Tdata));ddx_w0=zeros(N_dof,length(Tdata));
% %     k=ceil(j/2);
% %     for i=1:N_harm   % i=1,3,5
% %         x_w0(k,:)=x_w0(k,:)-index_global(i,j)*Harm_parameter_a(i,2*k-1)*Tdata.*sin(vector_w(i,k)*Tdata)+index_global(i,j)*Harm_parameter_a(i,2*k)*Tdata.*cos(vector_w(i,k)*Tdata);
% %         dx_w0(k,:)=dx_w0(k,:)-(index_global(i,j)*Harm_parameter_a(i,2*k-1)*sin(vector_w(i,k)*Tdata)+vector_w(i,k)*index_global(i,j)*Tdata*Harm_parameter_a(i,2*k-1).*cos(vector_w(i,k)*Tdata))+...
% %             (index_global(i,j)*Harm_parameter_a(i,2*k)*cos(vector_w(i,k)*Tdata)-vector_w(i,k)*index_global(i,j)*Tdata*Harm_parameter_a(i,2*k).*sin(vector_w(i,k)*Tdata));
% %         ddx_w0(k,:)=ddx_w0(k,:)-(2*vector_w(i,k)*index_global(i,j)*Harm_parameter_a(i,2*k-1)*cos(vector_w(i,k)*Tdata)-vector_w(i,k)^2*index_global(i,j)*Tdata*Harm_parameter_a(i,2*k-1).*sin(vector_w(i,k)*Tdata))-...
% %             (2*vector_w(i,k)*index_global(i,j)*Harm_parameter_a(i,2*k)*sin(vector_w(i,k)*Tdata)+vector_w(i,k)^2*index_global(i,j)*Tdata*Harm_parameter_a(i,2*k).*cos(vector_w(i,k)*Tdata));
% %     end
% %     k=ceil(j/2);
%     for ij=1:N_dof
%         for i=1:N_harm   % i=1,3,5
%             x_w0(ij,:)=x_w0(ij,:)-index_global(i,j)*Harm_parameter_a(i,2*ij-1)*Tdata.*sin(vector_w(i,ij)*Tdata)+index_global(i,j)*Harm_parameter_a(i,2*ij)*Tdata.*cos(vector_w(i,ij)*Tdata);
%             dx_w0(ij,:)=dx_w0(ij,:)-(index_global(i,j)*Harm_parameter_a(i,2*ij-1)*sin(vector_w(i,ij)*Tdata)+vector_w(i,ij)*index_global(i,j)*Tdata*Harm_parameter_a(i,2*ij-1).*cos(vector_w(i,ij)*Tdata))+...
%                 (index_global(i,j)*Harm_parameter_a(i,2*ij)*cos(vector_w(i,ij)*Tdata)-vector_w(i,ij)*index_global(i,j)*Tdata*Harm_parameter_a(i,2*ij).*sin(vector_w(i,ij)*Tdata));
%             ddx_w0(ij,:)=ddx_w0(ij,:)-(2*vector_w(i,ij)*index_global(i,j)*Harm_parameter_a(i,2*ij-1)*cos(vector_w(i,ij)*Tdata)-vector_w(i,ij)^2*index_global(i,j)*Tdata*Harm_parameter_a(i,2*ij-1).*sin(vector_w(i,ij)*Tdata))-...
%                 (2*vector_w(i,ij)*index_global(i,j)*Harm_parameter_a(i,2*ij)*sin(vector_w(i,ij)*Tdata)+vector_w(i,ij)^2*index_global(i,j)*Tdata*Harm_parameter_a(i,2*ij).*cos(vector_w(i,ij)*Tdata));
%         end
%     end
%     residual(N_dof*j+1:(1+j)*N_dof,:)=M*ddx_w0+C*dx_w0+K*x_w0+3*K3q*(x.^2.*x_w0);
% end
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
    residual(N_dof*(i)+1:N_dof*(i+1),:)=ddx_a+e*dx_a-e*(2*x.*dx.*x_a+dx_a.*x.^2)+w0_2*x_a+3*a*(x.^2.*x_a);
end

residual=residual';

