clc
clear all
close all
format long
syms q;
s=tf('s');
%% system 1 disrete
G1=2/((s^3)*(s+1));
% sampling time
Ts=0.04;
G1d=c2d(G1,Ts);
%% define idpoly system
%G num & den
A=G1d.den{1,1};
B=G1d.num{1,1};
%% input definning
U1=ones(1,250);
U2=-1*ones(1,250);
Uc=[U1 U2 U1 U2];



%% sys identify Projection algorithm
y=zeros(1,length(Uc));
yg=zeros(1,length(Uc));
Y=zeros(1,length(Uc));
% U=zeros(1,length(Uc));
 phi= zeros(1,8);
 alpha =1e6;
 p(:,:,4) = alpha*eye(length(phi(1,:)));
 teta =repmat([-3.95 5.85 -3.85 0.95 0.0000002 0.00002 0.000002 0.0000002 ]',1,4);
%  teta=[zeros(1,4);zeros(1,4)]';
% phi= zeros(1,);
% teta = [1 1 1 1;1 1 1 1]';
% K = zeros(length(phi_H(1,:)));
% P = zeros(length(phi_H(1,:)));
% alpha =1e10;
% P(:,:,3) = alpha*eye(length(phi_H(1,:)));
% teta_H_1 = [zeros(6,1),zeros(6,1),zeros(6,1)];
A_hat=zeros(1,4);
B_hat=zeros(1,4);
% gama=1.25;
fi=0.94;
c1=1;
c2=0.001;
c3=0.001;
e=zeros(1,length(Uc));
a=zeros(1,length(Uc));
U=zeros(1,length(Uc));


for i=5:length(Uc)
    if i<=15
 y(i)=-(A(2)*y(i-1)+A(3)*y(i-2)+A(4)*y(i-3))+B(2)*Uc(i-1)+B(3)*Uc(i-2)+B(4)*Uc(i-3)+B(5)*Uc(i-4);
 phi(i,:)=[-y(i-1) -y(i-2) -y(i-3) -y(i-4) Uc(i-1) Uc(i-2)  Uc(i-3)  Uc(i-4) ];
y_e(i)=phi(i-1,:)*teta(:,i-1);
e(i)=y(i)-y_e(i);
etha_plus(i)=c1*max(Uc(i))+c2+c3;
g(i)=etha_plus(i)/((1-fi)^(1/2));
if e(i)>g(i)
    f(i)=e(i)-g(i);
elseif e(i)<g(i) && e(i)>-g(i)
    f(i)=0;
elseif e(i)<-g(i)
    f(i)=e(i)+g(i);
end
if e(i)<=g(i) && e(i)>=-g(i)
    a(i)=0;
else
   a(i)=(fi*f(i))/e(i);
end
p(:,:,i-1)=p(:,:,i-2)-(a(i)*(p(:,:,i-2)*phi(i-1,:)'*phi(i-1,:)*p(:,:,i-2)/(1+phi(i-1,:)*p(:,:,i-2)*phi(i-1,:)')));

teta(:,i)=teta(:,i-1)+a(i)*(p(:,:,i-2)*phi(i,:)'*e(i))/(1+phi(i,:)*p(:,:,i-2)*phi(i,:)');
A_hat(i,1:4)=teta(1:4,i);
B_hat(i,1:4)=teta(5:8,i);
    end
if i>15


 y(i)=-(A(2)*y(i-1)+A(3)*y(i-2)+A(4)*y(i-3))+B(2)*U(i-1)+B(3)*U(i-2)+B(4)*U(i-3)+B(5)*U(i-4);
 phi(i,:)=[-y(i-1) -y(i-2) -y(i-3) -y(i-4) U(i-1) U(i-2)  U(i-3)  U(i-4) ];
y_e(i)=phi(i-1,:)*teta(:,i-1);
e(i)=y(i)-y_e(i);
etha_plus(i)=c1*max(Uc(i))+c2+c3;
g(i)=etha_plus(i)/((1-fi)^(1/2));
if e(i)>g(i)
    f(i)=e(i)-g(i);
elseif e(i)<g(i) && e(i)>-g(i)
    f(i)=0;
elseif e(i)<-g(i)
    f(i)=e(i)+g(i);
end
if e(i)<=g(i) && e(i)>=-g(i)
    a(i)=0;
else
   a(i)=(fi*f(i))/e(i);
end
p(:,:,i-1)=p(:,:,i-2)-(a(i)*(p(:,:,i-2)*phi(i-1,:)'*phi(i-1,:)*p(:,:,i-2)/(1+phi(i-1,:)*p(:,:,i-2)*phi(i-1,:)')));

teta(:,i)=teta(:,i-1)+a(i)*(p(:,:,i-2)*phi(i,:)'*e(i))/(1+phi(i,:)*p(:,:,i-2)*phi(i,:)');
A_hat(i,1:4)=teta(1:4,i);
B_hat(i,1:4)=teta(5:8,i);
end
%% closed loop
a1=0.0001;a2=0.002;
Ao_q=(q+a1)*(q+a2);
A_star=((1-0.9*q^-1)^5)*q^5*Ao_q;
A_Star=sym2poly(A_star);
T=sym2poly(Ao_q);



%% silvester matrix and Diopgh equation

E=[A_hat(i,4) 0 0 0 B_hat(i,4) 0 0 0;
    A_hat(i,3) A_hat(i,4) 0 0 B_hat(i,3) B_hat(i,4) 0 0;
    A_hat(i,2) A_hat(i,3) A_hat(i,4) 0 B_hat(i,2) B_hat(i,3) B_hat(i,4) 0;
    A_hat(i,1) A_hat(i,2) A_hat(i,3) A_hat(i,4) B_hat(i,1) B_hat(i,2) B_hat(i,3) B_hat(i,4);
    1 A_hat(i,1) A_hat(i,2) A_hat(i,3) 0 B_hat(i,1) B_hat(i,2) B_hat(i,3);
    0 1 A_hat(i,1) A_hat(i,2) 0 0 B_hat(i,1) B_hat(i,2);
    0 0 1 A_hat(i,1) 0 0 0 B_hat(i,1);
    0 0 0  1 0 0 0 0;
    ];

%% controller parameter
LP=E\A_Star';
L=LP(1:4,1)';
Lp=1-L;
Lp=flipdim(Lp,2);
P=LP(5:8,1)';

Y(i)=-(A_hat(i,1)*Y(i-1)+A_hat(i,2)*Y(i-2)+A_hat(i,3)*Y(i-3)+A_hat(i,4)*Y(i-4))+B_hat(i,1)*U(i-1)+B_hat(i,2)*U(i-2)+B_hat(i,3)*U(i-3)+B_hat(i,4)*U(i-4);

U(i)=-(L(2)*U(i-1)+L(3)*U(i-2)+L(4)*U(i-3))+P(1)*Uc(i)+P(2)*Uc(i-1)+P(3)*Uc(i-2)+P(4)*Uc(i-3)-P(1)*Y(i)-P(2)*Y(i-1)-P(3)*Y(i-2)-P(4)*Y(i-3);

if U(i)>2
    U(i)=2;
elseif U(i)<-2
    U(i)=-2;
end

end
%% silvester matrix and Diopgh equation

E_n=[ A(5)       0              0         0           B(5)        0          0           0;
    A(4)      A(5)            0         0           B(4)        B(5)       0           0;
    A(3)      A(4)            A(5)      0           B(3)        B(4)       B(5)        0;
    A(2)      A(3)            A(4)      A(5)        B(2)        B(3)       B(4)        B(5);
    1         A(2)            A(3)      A(4)        0           B(2)       B(3)        B(4);
    0         1               A(2)      A(3)        0           0          B(2)        B(3);
    0         0               1         A(2)        0           0          0           B(2);
    0         0               0         1           0           0          0           0;
    ];

%% controller parameter
LP_n=E_n\A_Star';
L_n=LP_n(1:4,1)';
Lp_n=1-L_n;
Lp_n=flipdim(Lp_n,2);
P_n=LP_n(5:8,1)';

Ug=zeros(1,length(Uc));
yg=zeros(1,length(Uc));
for i=5:length(Uc)
     yg(i)=-(A(2)*yg(i-1)+A(3)*yg(i-2)+A(4)*yg(i-3))+B(2)*Ug(i-1)+B(3)*Ug(i-2)+B(4)*Ug(i-3)+B(5)*Ug(i-4);

    
Ug(i)=(Lp_n(2)*Ug(i-1)+Lp_n(3)*Ug(i-2)+Lp_n(4)*Ug(i-3))+P_n(1)*Uc(i)+P_n(2)*Uc(i-1)+P_n(3)*Uc(i-2)+P_n(4)*Uc(i-3)-P_n(1)*yg(i)-P_n(2)*yg(i-1)-P_n(3)*yg(i-2)-P_n(4)*yg(i-3);

if Ug(i)>2
    Ug(i)=2;
elseif Ug(i)<-2
    Ug(i)=-2;
end
end

%% plotting
 figure;
 subplot(2,1,1)
 plot(Uc,'b','linewidth',3);hold on
 plot(Y,'g-','linewidth',3);hold on
 plot(0.95*Y,'r--','linewidth',2);
 legend('Uc','y1','y2');grid on
% hold on
%  plot(2.5*yg,'r--','linewidth',2.5);legend('Uc','Y','Yd');axis([0 500 -4 6]);grid on
subplot(2,1,2);title('control signal')
plot(0.9*U,'b','linewidth',2.5);hold on
plot(U,'r--','linewidth',1.5)
legend('u1','u2')
grid on