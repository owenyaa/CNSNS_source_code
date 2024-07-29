function yprime=ode_floquet(t1,y)
global A parameter_a M N_dof w N_harm ep k_n f1
Harm_parameter_a=parameter_a;
x=zeros(2,1);
for j=1:N_dof
    for i=1:N_harm   % i=1,3,5
        x(j,:)=x(j,:)+Harm_parameter_a(i,2*j-1)*cos((2*i-1)*w*t1)+Harm_parameter_a(i,2*j)*sin((2*i-1)*w*t1);
    end
end
df=[-ep*f1*cos(w*t1)+3*ep*k_n*(x(1,:)-x(2,:)).^2;3*ep*k_n*(x(2,:)-x(1,:)).^2];
% df=[w*ep*f1*sin(w*t1)+3*ep*k_n*(x(1,:)-x(2,:)).^2;3*ep*k_n*(x(2,:)-x(1,:)).^2];
% df=[3*ep*k_n*(x(1,:)-x(2,:)).^2;3*ep*k_n*(x(2,:)-x(1,:)).^2];
df=[0;0;-M\df];
temp=zeros(4,4);temp(:,2)=df;


% %.....q'=Aq+V*K3q(3,3)*beta^3
% %.....dq'=A*dq+V*3*K3q(3,3)*beta^2
% 
% temp=[1;1];%.......Mq''+Cq'+Kq+temp*K3q(3,3)*beta^3=0
% V=-[0;0;M\temp];
% temp=zeros(4,4);
% % temp(:,2)=V*3*K3q(3,3)*beta_t^2;
% temp(:,2)=[0;0;-ep*f1*cos(w*t1)+3*ep*k_n*(x(1,:)-x(2,:)).^2;3*ep*k_n*(x(2,:)-x(1,:)).^2].*V;

yprime=(A+temp)*y;