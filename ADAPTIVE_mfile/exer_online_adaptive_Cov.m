%% Exer_2_adaptive_online_H_parameter changing slowly_Cov resetiing
clc
clear all
close all
format long
%% define tranfer function
s=tf('s');
G=(s+0.5)/((3*s+1)*(s^2+0.7*s+1));
H=(s+0.5)/((3*s-1)*(s^2+0.7*s+1));
%% discrete system with 20 times bandwidth sapmilg frequency
fg=bandwidth(G);
fh=bandwidth(H);
Tg=1/(20*fg);
Th=1/(20*fh);
Gd=c2d(G,Tg);
Hd=c2d(H,Th);
%% define idpoly system
%G num & den
Ag=Gd.den{1,1};
Bg=Gd.num{1,1};
AG=[Ag(1,2:4) Bg(1,2:4)]';
% idpoly G system
m0_G=idpoly(Ag,Bg)
% H num & den
Ah=Hd.den{1,1};
Bh=Hd.num{1,1};
AH=[Ah(1,2:4) Bh(1,2:4)]';
% idpoly H system
m0_H=idpoly(Ah,Bh)

%% system H with random inputs
n=1000;
u_1=[ zeros(1,3) rand(1,n)]';
yh=zeros(n,1);
% u = [0;0;u_1];
phi_H= zeros(1,6);

k = zeros(length(phi_H(1,:)));
p = zeros(length(phi_H(1,:)));
alpha =1e10;
p(:,:,3) = alpha*eye(length(phi_H(1,:)));
teta_H = [zeros(6,1),zeros(6,1),zeros(6,1)];
t=0.0001:-0.00000001:-0.000001*size(u_1);
t_1=0.01:-0.00001:-0.0001*size(u_1);

%% adaptive identification system H
for i=4:1:max(size(u_1))
           A_h(1,i)=Ag(1,2)+t_1(i);
    A_h(2,i)=Ag(1,3)+t_1(i);
    A_h(3,i)=Ag(1,4)+t_1(i);
    B_h(1,i)=Bg(1,2)+t(i);
    B_h(2,i)=Bg(1,3)+t(i);
    B_h(3,i)=Bg(1,4)+t(i);
 yh(i)=-((Ah(1,2)+t_1(i))*yh(i-1))-(Ah(1,3)+t_1(i))*yh(i-2)-((Ah(1,4)+t_1(i))*yh(i-3))+((Bh(1,2)+t(i))*u_1(i-1))+((Bh(1,3)+t(i))*u_1(i-2))+((Bh(1,4)+t(i))*u_1(i-3));   
 phi_H(i,:)=[-yh(i-1) -yh(i-2) -yh(i-3) u_1(i-1) u_1(i-2) u_1(i-3)];
  if (i/20)==round(i/20);
      k=zeros(length(phi_H(i,:)));
    p = zeros(length(phi_H(i,:)));
alpha =1e14;
p(:,:,i) = alpha*eye(length(phi_H(i,:)));
  end
    k(:,i)=p(:,:,i-1)*phi_H(i,:)'/((1)+phi_H(i,:)*p(:,:,i-1)*phi_H(i,:)');
    p(:,:,i)=(eye(length(phi_H(i,:)))-k(:,i)*phi_H(i,:))*p(:,:,i-1);
    teta_H(:,i)=(teta_H(:,i-1)+k(:,i)*(yh(i)-phi_H(i,:)*teta_H(:,i-1)));

 
end
init ='z';
yu = sim (m0_H,u_1);
%% plotting system parameters
figure;

subplot(2,3,1)
plot(4:i,A_h(1,i),'r*');
hold on
plot(teta_H(1,:),'b--');legend('a1')
subplot(2,3,2)
plot(4:i,A_h(2,i),'r*')
hold on
plot(teta_H(2,:),'b--');legend('a2')
subplot(2,3,3)
plot(4:i,A_h(3,i),'r*')
hold on
plot(teta_H(3,:),'b--');legend('a3')
subplot(2,3,4)
plot(4:i,B_h(1,i),'r*')
hold on
plot(teta_H(4,:),'b--');legend('b1')
subplot(2,3,5)
plot(4:i,B_h(2,i),'r*')
hold on
plot(teta_H(5,:),'b--');legend('b2')
subplot(2,3,6)
plot(4:i,B_h(3,i),'r*')
hold on
plot(teta_H(6,:),'b--');legend('b3')
%  figure
%  plot(u_1,'r');legend('random input')
figure
subplot(2,1,1)
plot(yu,'b');legend('real output')
subplot(2,1,2)
plot(phi_H*teta_H,'black--')
;legend('phi*teta')