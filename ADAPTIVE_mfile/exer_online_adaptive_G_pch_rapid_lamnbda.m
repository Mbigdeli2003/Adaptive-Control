%% exer_online_adaptive_G_pch_rapid_lambda
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
%% system G with random inputs
n=1000;
u_1=[ zeros(1,3) rand(1,n)]';
yg=zeros(n,1);
% u = [0;0;u_1];
phi_G= zeros(1,6);
k = zeros(length(phi_G(1,:)));
p = zeros(length(phi_G(1,:)));
alpha =1e10;
p(:,:,3) = alpha*eye(length(phi_G(1,:)));
teta_G = [zeros(6,1),zeros(6,1),zeros(6,1)];

Betao=0.05*rand(1,n+3);
alfa=10*rand(1,n+3);
ch=0.0002*ones(n+3,1);
ch_1=0.0000002*ones(n+3,1);


%% adaptive identification system G
for i=4:1:max(size(u_1))
    lambda(i)=1-(Betao(i)*exp(-alfa(i)*i));
 yg(i)=-((Ag(1,2))*yg(i-1))-(Ag(1,3))*yg(i-2)-((Ag(1,4))*yg(i-3))+((Bg(1,2))*u_1(i-1))+((Bg(1,3))*u_1(i-2))+((Bg(1,4))*u_1(i-3));   
if(i/100)>=1
    
 yg(i)=-((Ag(1,2)+ch(i))*yg(i-1))-(Ag(1,3)+ch(i))*yg(i-2)-((Ag(1,4)+ch(i))*yg(i-3))+((Bg(1,2)+ch_1(i))*u_1(i-1))+((Bg(1,3)+ch_1(i))*u_1(i-2))+((Bg(1,4)+ch_1(i))*u_1(i-3));    
 
  A_g(1,i)=Ag(1,2)+ch(i);
    A_g(2,i)=Ag(1,3)+ch(i);
    A_g(3,i)=Ag(1,4)+ch(i);
    B_g(1,i)=Bg(1,2)+ch_1(i);
    B_g(2,i)=Bg(1,3)+ch_1(i);
    B_g(3,i)=Bg(1,4)+ch_1(i);
end
 phi_G(i,:)=[-yg(i-1) -yg(i-2) -yg(i-3) u_1(i-1) u_1(i-2) u_1(i-3)];
    k(:,i)=p(:,:,i-1)*phi_G(i,:)'/((lambda(i)*1)+phi_G(i,:)*p(:,:,i-1)*phi_G(i,:)');
    p(:,:,i)=(eye(length(phi_G(i,:)))-k(:,i)*phi_G(i,:))*p(:,:,i-1);
    teta_G(:,i)=(teta_G(:,i-1)+k(:,i)*(yg(i)-phi_G(i,:)*teta_G(:,i-1)))/lambda(i);

end
init ='z';
yu = sim (m0_G,u_1);
%% plotting
figure;
subplot(2,3,1)
plot(4:i+3,A_g(1,:),'r*');
hold on
plot(teta_G(1,:),'b--');legend('a1')
subplot(2,3,2)
plot(4:i+3,A_g(2,:),'r*')
hold on
plot(teta_G(2,:),'b--');legend('a2')
subplot(2,3,3)
plot(4:i+3,A_g(3,:),'r*')
hold on
plot(teta_G(3,:),'b--');legend('a3')
subplot(2,3,4)
plot(4:i+3,B_g(1,:),'r*')
hold on
plot(teta_G(4,:),'b--');legend('b1')
subplot(2,3,5)
plot(4:i+3,B_g(2,:),'r*')
hold on
plot(teta_G(5,:),'b--');legend('b2')
subplot(2,3,6)
plot(4:i+3,B_g(3,:),'r*')
hold on
plot(teta_G(6,:),'b--');legend('b3')
%  figure
%  plot(u_1,'r');legend('random input')
figure
subplot(2,1,1)
plot(yu,'b');legend('real output')
subplot(2,1,2)
plot(phi_G*teta_G,'black--')
;legend('phi*teta')





    
    
    
    
    















