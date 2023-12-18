%% STR Direct_nozero pole cancelled-martabe3-Delay
clc
clear all
close all
syms q
%% deifine Uc
%define T 3Ts=3*8.93
t=0:150;
% Uc=zeros(9*length(T));
u1=2*ones(1,length(t));
u2=-1*ones(1,length(t));
u3=ones(1,length(t));
Uc=[zeros(1,3) u1 u2 u3];


%% define Tf
s=tf('s');
%unstable system
H=(s+0.5)/((3*s-1)*(s^2+1*s+1.2));
%desired system
G=(s+0.5)/((3*s+1)*(s^2+1*s+1.2));
%% dicrtisize sys with 20*bandwidth
fh=bandwidth(H);
fg=bandwidth(G);
Th=0.2;
Hd=c2d(H,Th);
Tg=Th;
Gd=c2d(G,Tg);
% H num & den
Ah=Hd.den{1,1};
Bh=Hd.num{1,1};
B=Bh(2:end);
Ag=Gd.den{1,1};
Am=Ag;
Bg=Gd.num{1,1};
B_g=Bg(2:end);
Bm=B_g;
%% Direct STR
Yf=[ zeros(1,length(Uc))];
yg=[ zeros(1,length(Uc))];
Y=[ zeros(1,length(Uc))];
Uf=[ zeros(1,length(Uc))];
U=[ zeros(1,length(Uc))];
yh=[ zeros(1,length(Uc))];
a1=0.01;a2=0.02;
A0_q=(q+a1)*(q+a2);
Ao=sym2poly(A0_q);
Am_q=poly2sym(Am,q);
Ac_q=A0_q*Am_q;
Ac=sym2poly(Ac_q);


%% system identifying
phi_H= zeros(1,6);

k = zeros(length(phi_H(1,:)));
p = zeros(length(phi_H(1,:)));
alpha =1e10;
p(:,:,5) = alpha*eye(length(phi_H(1,:)));
teta_H = [zeros(6,1),zeros(6,1),zeros(6,1),zeros(6,1),zeros(6,1)];
phi_H= zeros(1,6);

K = zeros(length(phi_H(1,:)));
P = zeros(length(phi_H(1,:)));
alpha =1e10;
P(:,:,5) = alpha*eye(length(phi_H(1,:)));
teta_H_1 = [zeros(6,1),zeros(6,1),zeros(6,1),zeros(6,1),zeros(6,1)];
A=zeros(1,3);
B=zeros(1,3);






%% controller identifying

phi=0.01* ones(1,6);

k = zeros(length(phi_H(1,:)));
p = zeros(length(phi_H(1,:)));
alpha =1e20;
p(:,:,5) = alpha*eye(length(phi(1,:)));
% teta = [zeros(6,1),zeros(6,1),zeros(6,1)];
phi=0.01* ones(1,6);

K = zeros(length(phi(1,:)));
P = zeros(length(phi(1,:)));
alpha =1e20;
P(:,:,5) = alpha*eye(length(phi(1,:)));
teta = [0.01*ones(6,1),0.01*ones(6,1),0.01*ones(6,1),0.01*ones(6,1),0.01*ones(6,1)];
R=zeros(length(Uc),3);
S=zeros(length(Uc),3);



%% adaptive identification system H
% for o=4:length(Uc)

for i=6:length(Uc)
    %desired system
    
    yg(i)=-(Ag(1,2)*yg(i-1))-Ag(1,3)*yg(i-2)-(Ag(1,4)*yg(i-3))+(Bg(1,2)*Uc(i-1))+(Bg(1,3)*Uc(i-2))+(Bg(1,4)*Uc(i-3));
  if i<30
    Y(i)=-(Ah(1,2)*Y(i-3))-(Ah(1,3)*Y(i-4))-(Ah(1,4)*Y(i-5))+(Bh(1,2)*Uc(i-3))+(Bh(1,3)*Uc(i-4))+(Bh(1,4)*Uc(i-5));
  phi_H(i,:)=[-Y(i-3) -Y(i-4) -Y(i-5) Uc(i-3) Uc(i-4) Uc(i-5)];
  elseif i>=30
    
       Y(i)=-(A(i-1,1)*Y(i-3))-(A(i-1,2)*Y(i-4))-(A(i-1,3)*Y(i-5))+(B(i-1,1)*U(i-3))+(B(i-1,2)*U(i-4))+(B(i-1,3)*U(i-5)); 
 phi_H(i,:)=[-Y(i-3) -Y(i-4) -Y(i-5) U(i-3) U(i-4) U(i-5)];
  end
  
    k(:,i)=p(:,:,i-1)*phi_H(i,:)'/(1+phi_H(i,:)*p(:,:,i-1)*phi_H(i,:)');
    p(:,:,i)=(eye(length(phi_H(i,:)))-k(:,i)*phi_H(i,:))*p(:,:,i-1);
    teta_H(:,i)=teta_H(:,i-1)+k(:,i)*(Y(i)-phi_H(i,:)*teta_H(:,i-1));
A(i,1:3)=teta_H(1:3,i);
B(i,1:3)=teta_H(4:6,i);
    
    
   
        phi(i,:)=[Uf(i-3) Uf(i-4) Uf(i-5) Yf(i-3) Yf(i-4) Yf(i-5)];
    k(:,i)=p(:,:,i-1)*phi(i,:)'/(1+phi(i,:)*p(:,:,i-1)*phi(i,:)');
    p(:,:,i)=(eye(length(phi_H(i,:)))-k(:,i)*phi_H(i,:))*p(:,:,i-1);
    teta(:,i)=teta(:,i-1)+k(:,i)*(Y(i)-phi(i,:)*teta(:,i-1));
R(i,1:3)=teta(1:3,i);
S(i,1:3)=teta(4:6,i);
t0=0.401*(sum(Am))/(sum(B(i,1:3)));
T=t0*Ao;
 Uf(i)=B(i,1)*U(i-1)+B(i,2)*U(i-2)+B(i,3)*U(i-3)-(Ac(1,2)*Uf(i-1)+Ac(1,3)*Uf(i-2)+Ac(1,4)*Uf(i-3)+Ac(1,5)*Uf(i-4)+Ac(1,6)*Uf(i-5));
 Yf(i)=B(i,1)*Y(i-1)+B(i,2)*Y(i-2)+B(i,3)*Y(i-3)-(Ac(1,2)*Yf(i-1)+Ac(1,3)*Yf(i-2)+Ac(1,4)*Yf(i-3)+Ac(1,5)*Yf(i-4)+Ac(1,6)*Yf(i-5));
U(i)=((T(1,1)*Uc(i)+T(1,2)*Uc(i-1)+T(1,3)*Uc(i-2))-(S(i,1)*Y(i)+S(i,2)*Y(i-1)+S(i,3)*Y(i-2))-(R(i,2)*U(i-1)+R(i,3)*U(i-2)))/(R(i,1));
    
%  Y(i)=R(1,1)*Uf(i-1)+R(1,2)*Uf(i-2)+R(1,3)*Uf(i-3)+S(1,1)*Uf(i-1)+S(1,2)*Uf(i-2)+S(1,3)*Uf(i-3);
end
%% plotting
 figure;
 subplot(2,1,1)
 plot(Uc,'b','linewidth',3);axis([0 500 -4 6]);hold on
 plot(2.5*Y,'g-','linewidth',3);axis([0 500 -4 6]);hold on
 plot(2.5*yg,'r--','linewidth',2.5);legend('Uc','Y','Yd');axis([0 500 -4 6]);grid on
subplot(2,1,2);
plot(U,'b','linewidth',1);legend('control signal');;axis([0 500 -4 6])
grid on
TT=repmat(T,456,1);
figure;
plot(R(1:i,2)./R(1:i,1),'b-^');axis([0 500 -1.4 1.4]);hold on
plot(R(1:i,3)./R(1:i,1),'b-*');;axis([0 500 -1.4 1.4]);hold on
plot(TT(1:i,1)./R(1:i,1),'b-.');axis([0 500 -1.4 1.4]);hold on
plot(TT(1:i,2)./R(1:i,1),'g--');axis([0 500 -1.4 1.4]);hold on
plot(TT(1:i,3)./R(1:i,1),'r-');axis([0 500 -1.4 1.4]);legend('r1/r0','r2/r0','t0/r0','t1/r0','t2/r0','t3/r0')


