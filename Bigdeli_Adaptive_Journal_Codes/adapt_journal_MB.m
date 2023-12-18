%%example one
clc
clear all
close all
format short
syms q;
s=tf('s');
%% system 1 disrete
G1=5/(s*(s+1));
% sampling time
Ts=0.04;
G1d=c2d(G1,Ts);
%% define idpoly system
%G num & den
A1=G1d.den{1,1};
B1=G1d.num{1,1};
%% input definning
U1=ones(1,250);
U2=-1*ones(1,250);
Uc=[U1 U2 U1 U2];

%% sys identify Projection algorithm
y=zeros(1,length(Uc));
yg=zeros(1,length(Uc));
Y=zeros(1,length(Uc));
% U=zeros(1,length(Uc));
 phi= zeros(1,4);
 alpha =1e10;
 p(:,:,2) = alpha*eye(length(phi(1,:)));
 teta = [-1.9 0.8 0.001 0.001;-1.9 0.9 0.001 0.001]';
%  teta=[zeros(1,4);zeros(1,4)]';
phi= zeros(1,4);
% teta = [1 1 1 1;1 1 1 1]';
% K = zeros(length(phi_H(1,:)));
% P = zeros(length(phi_H(1,:)));
% alpha =1e10;
% P(:,:,3) = alpha*eye(length(phi_H(1,:)));
% teta_H_1 = [zeros(6,1),zeros(6,1),zeros(6,1)];
A_hat=zeros(1,2);
B_hat=zeros(1,2);
% gama=1.25;
fi=0.7;
c1=1;
c2=0.001;
c3=0.001;
e=zeros(1,length(Uc));
a=zeros(1,length(Uc));
U=zeros(1,length(Uc));
for i=3:length(Uc)
    if i<=10
 y(i)=-(A1(2)*y(i-1)+A1(3)*y(i-2))+B1(2)*Uc(i-1)+B1(3)*Uc(i-2);
 phi(i,:)=[-y(i-1) -y(i-2) Uc(i-1) Uc(i-2)];
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
A_hat(i,1:2)=teta(1:2,i);
B_hat(i,1:2)=teta(3:4,i);
    end
if i>10
     y(i)=-(A1(2)*y(i-1)+A1(3)*y(i-2))+B1(2)*U(i-1)+B1(3)*U(i-2);
 phi(i,:)=[-y(i-1) -y(i-2) U(i-1) U(i-2)];
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
A_hat(i,1:2)=teta(1:2,i);
B_hat(i,1:2)=teta(3:4,i);
end   
    
%% closed loop
a=0.001;
Ao_q=q+a;
A_star=((1-0.9*q^-1)^2)*q^2*Ao_q;
A_Star=sym2poly(A_star);
T=sym2poly(Ao_q);

%% slivester matrix & Dioph equation
E=[ A_hat(i,2)     0           B_hat(i,2)      0 ;    
    A_hat(i,1)    A_hat(i,2)  B_hat(i,1)      B_hat(i,2);
    1             A_hat(i,1)  0               B_hat(i,1);
    0             1           0               0            ];


LP=E\A_Star';
L=LP(1:2,1)';
Lp=1-L;
Lp=flipdim(Lp,2);
P=LP(3:4,1)';
Y(i)=-(A_hat(i,1)*Y(i-1)+A_hat(i,2)*Y(i-2))+B_hat(i,1)*U(i-1)+B_hat(i,2)*U(i-2);

U(i)=Lp(2)*U(i-1)+T(1)*Uc(i)+T(2)*Uc(i-1)-P(1)*Y(i)-P(2)*Y(i-1);

if U(i)>1
    U(i)=1;
elseif U(i)<-1
    U(i)=-1;
end


end
%% non adptive case
% closed loop
a=0.001;
Ao_q=q+a;
A_star=((1-0.9*q^-1)^2)*q^2*Ao_q;
A_Star=sym2poly(A_star);
T=sym2poly(Ao_q);

%% slivester matrix & Dioph equation
E_n=[ A1(3)     0           B1(3)    0 ;    
    A1(2)       A1(3)       B1(2)    B1(3);
    1           A1(2)       0        B1(2);
    0           1           0        0            ];


LP_n=E_n\A_Star';
L_n=LP_n(1:2,1)';
Lp_n=1-L_n;
Lp_n=flipdim(Lp_n,2);
P_n=LP_n(3:4,1)';

Ug=zeros(1,length(Uc));
yg=zeros(1,length(Uc));
for i=3:length(Uc)
 yg(i)=-(A1(2)*yg(i-1)+A1(3)*yg(i-2))+B1(2)*Ug(i-1)+B1(3)*Ug(i-2);
 
Ug(i)=Lp_n(2)*Ug(i-1)+T(1)*Uc(i)+T(2)*Uc(i-1)-P_n(1)*yg(i)-P_n(2)*yg(i-1);

if Ug(i)>1
    Ug(i)=1;
elseif Ug(i)<-1
    Ug(i)=-1;
end



end

%% plotting
 figure;
 subplot(2,1,1)
 plot(Uc,'b','linewidth',3);hold on
 plot(5*Y,'g-','linewidth',3);hold on
 plot(1.3*yg,'r--','linewidth',2);
 legend('Uc','y1','y2');grid on
% hold on
%  plot(2.5*yg,'r--','linewidth',2.5);legend('Uc','Y','Yd');axis([0 500 -4 6]);grid on
subplot(2,1,2);title('control signal')
plot(U,'b','linewidth',2.5);hold on
plot(Ug,'r--','linewidth',1.5)
legend('u1','u2')
grid on







