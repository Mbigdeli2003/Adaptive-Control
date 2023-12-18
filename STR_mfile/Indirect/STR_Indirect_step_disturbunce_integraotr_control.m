%% STR-dynamic-Indirect feedback-pole placemnet-nozero cancelled_white
%% step disturbunce integretor control
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
Bg=Gd.num{1,1};
B_g=Bg(2:end);
%system polynaminals defining with q

%% system identifying
%% sys identify
% noise loading, minor(var=0.01),major(var=2),average(var=0.1)
load noise_major
load noise_av
load noise_minor
yh=zeros(1,length(Uc));
yg=zeros(1,length(Uc));
Y=zeros(1,length(Uc));
U=zeros(1,length(Uc));
phi_H= zeros(1,6);

k = zeros(length(phi_H(1,:)));
p = zeros(length(phi_H(1,:)));
alpha =1e10;
p(:,:,3) = alpha*eye(length(phi_H(1,:)));
teta_H = [zeros(6,1),zeros(6,1),zeros(6,1)];
phi_H= zeros(1,6);

K = zeros(length(phi_H(1,:)));
P = zeros(length(phi_H(1,:)));
alpha =1e10;
P(:,:,3) = alpha*eye(length(phi_H(1,:)));
teta_H_1 = [zeros(6,1),zeros(6,1),zeros(6,1)];
A=zeros(1,3);
B=zeros(1,3);
%% defining step disturbunce
e=0.8*ones(1,2*length(Uc));
t=0;
%% adaptive identification system H
% for o=4:length(Uc)

for i=4:length(Uc)
    t=t+1;
    
%      if  teta_H(1:3,i-1)'-Ah(1:3)==0 & teta_H(4:6,i-1)'-Bh(1:3)==0
 yh(i)=-(Ah(1,2)*yh(i-1))-Ah(1,3)*yh(i-2)-(Ah(1,4)*yh(i-3))+(Bh(1,2)*Uc(i-1))+(Bh(1,3)*Uc(i-2))+(Bh(1,4)*Uc(i-3)); 
 if t>=160
  yh(i)=-(Ah(1,2)*yh(i-1))-Ah(1,3)*yh(i-2)-(Ah(1,4)*yh(i-3))+(Bh(1,2)*Uc(i-1))+(Bh(1,3)*Uc(i-2))+(Bh(1,4)*Uc(i-3))+e(t);
 end
  phi_H(i,:)=[-yh(i-1) -yh(i-2) -yh(i-3) Uc(i-1) Uc(i-2) Uc(i-3)];
    k(:,i)=p(:,:,i-1)*phi_H(i,:)'/(1+phi_H(i,:)*p(:,:,i-1)*phi_H(i,:)');
    p(:,:,i)=(eye(length(phi_H(i,:)))-k(:,i)*phi_H(i,:))*p(:,:,i-1);
    teta_H(:,i)=teta_H(:,i-1)+k(:,i)*(yh(i)-phi_H(i,:)*teta_H(:,i-1));
A(i,1:3)=teta_H(1:3,i);
B(i,1:3)=teta_H(4:6,i);
 
if i>=10
 %% define desired system and real system with pole zero cancellation 
 Am=Ag;
 Bm=B_g;
 Am_q=poly2sym(Am,q);
 %% control with diophantine  
 %observer
 a0=0.02;a1=0.01;%defining A0 coeffs
 Ao_q=(q+a0)*(q+a1);
 B_plus=1;  
 B_minus=B(i,1);
 B_plus_q=poly2sym(B_plus,q);
 B_q=poly2sym(B(i,1:3),q);
 A_q=poly2sym([1 A(i,1:3)],q);
   
   %% desired system poly
Ac_q=Ao_q*Am_q*B_plus_q;
Ac_c=coeffs(Ac_q,q);
Ac=sym2poly(Ac_c);

%% silvester matrix and Dioph equation
E=[A(i,3) 0 0 B(i,3) 0 0;...
   A(i,2) A(i,3) 0 B(i,2) B(i,3) 0;...
   A(i,1) A(i,2) A(i,3) B(i,1) B(i,2) B(i,3);...
   1 A(i,1) A(i,2) 0 B(i,1) B(i,2);...
   0 1 A(i,1) 0 0 B(i,1);...
   0 0 1 0 0 0];
% controller parameter
T=sym2poly(Ao_q);
RS=E\Ac';
R=RS(1:3,1)';
S=RS(4:6,1)';
R=flipdim(R,2);
S=flipdim(S,2);
% RI=[R 0];
% SI=[S 0];
%% integrator contrller step disturbunce
% if t>=160
X=[1 0.85];
X_q=poly2sym(X,q);
R_q=poly2sym(R,q);
S_q=poly2sym(S,q);
%new controller with integrator
R_Q=(X_q*R_q)+(sum(X)*(sum(R)/-sum(B(i,1:3)))*B_q);
S_Q=(X_q*S_q)-(sum(X)*(sum(R)/-sum(B(i,1:3)))*A_q);
RI=sym2poly(R_Q);
SI=sym2poly(S_Q);

%% testing controller obtaiend parameters Ru=TUc-Sy

    % main system
   yh(i)=-(Ah(1,2)*yh(i-1))-Ah(1,3)*yh(i-2)-(Ah(1,4)*yh(i-3))+(Bh(1,2)*Uc(i-1))+(Bh(1,3)*Uc(i-2))+(Bh(1,4)*Uc(i-3));

  
     
    Y(i)=-(A(i,1)*Y(i-1))-A(i,2)*Y(i-2)-(A(i,3)*Y(i-3))+(B(i,1)*U(i-1))+(B(i,2)*U(i-2))+(B(i,3)*U(i-3));
    
    U(i)=-(RI(1,2)*U(i-1)+RI(1,3)*U(i-2)+RI(1,4)*U(i-3))+(T(1,1)*Uc(i)+T(1,2)*Uc(i-1)+T(1,3)*Uc(i-2))-(SI(1,1)*Y(i)+SI(1,2)*Y(i-1)+SI(1,3)*Y(i-2)+SI(1,4)*Y(i-3));
  
  
  
  %  desired system
yg(i)=-(Ag(1,2)*yg(i-1))-Ag(1,3)*yg(i-2)-(Ag(1,4)*yg(i-3))+(Bg(1,2)*Uc(i-1))+(Bg(1,3)*Uc(i-2))+(Bg(1,4)*Uc(i-3));
end
end

% end
%% plotting
 figure;
 subplot(2,1,1)
 plot(Uc,'b','linewidth',3);axis([0 500 -4 6]);hold on
 plot(4*Y,'g-','linewidth',3);axis([0 500 -4 6]);hold on
 plot(2.5*yg,'r--','linewidth',2.5);legend('Uc','Y','Yd');axis([0 500 -4 6]);grid on
subplot(2,1,2);
plot(U,'b','linewidth',2.5);legend('control signal')
grid on
%% system ident plotting
figure;
subplot(2,3,1)
plot(4:i,Ah(1,2),'r*');
hold on
plot(teta_H(1,:),'b--');legend('a1')
subplot(2,3,2)
plot(4:i,Ah(1,3),'r*')
hold on
plot(teta_H(2,:),'b--');legend('a2')
subplot(2,3,3)
plot(4:i,Ah(1,4),'r*')
hold on
plot(teta_H(3,:),'b--');legend('a3')
subplot(2,3,4)
plot(4:i,Bh(1,2),'r*')
hold on
plot(teta_H(4,:),'b--');legend('b1')
subplot(2,3,5)
plot(4:i,Bh(1,3),'r*')
hold on
plot(teta_H(5,:),'b--');legend('b2')
subplot(2,3,6)
plot(4:i,Bh(1,4),'r*')
hold on
plot(teta_H(6,:),'b--');legend('b3')




















