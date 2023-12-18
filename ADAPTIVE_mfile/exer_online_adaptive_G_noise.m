%% Exer_2_adaptive_online_G_white noise
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
n=200;
u_1=[ zeros(1,3) ones(1,n)]';
varN=2;% noise variance
r=random('normal',0,varN,n,1);
noise=[zeros(1,3) r']';
c=[1 0.1 0.5 2 3];
yg=zeros(n,1);
% u = [0;0;u_1];
phi_G= zeros(1,6);
k = zeros(length(phi_G(1,:)));
p = zeros(length(phi_G(1,:)));
alpha =1e6;
p(:,:,3) = alpha*eye(length(phi_G(1,:)));
teta_G = [zeros(6,1),zeros(6,1),zeros(6,1)];


%% adaptive identification system G
for i=4:1:max(size(u_1))
 yg(i)=-(Ag(1,2)*yg(i-1))-Ag(1,3)*yg(i-2)-(Ag(1,4)*yg(i-3))+(Bg(1,2)*u_1(i-1))+(Bg(1,3)*u_1(i-2))+(Bg(1,4)*u_1(i-3))+c(1,1)*noise(i);   
  phi_G(i,:)=[-yg(i-1) -yg(i-2) -yg(i-3) u_1(i-1) u_1(i-2) u_1(i-3)];
    k(:,i)=p(:,:,i-1)*phi_G(i,:)'/(1+phi_G(i,:)*p(:,:,i-1)*phi_G(i,:)');
    p(:,:,i)=(eye(length(phi_G(i,:)))-k(:,i)*phi_G(i,:))*p(:,:,i-1);
    teta_G(:,i)=teta_G(:,i-1)+k(:,i)*(yg(i)-phi_G(i,:)*teta_G(:,i-1));
end
init ='z';
yu = sim (m0_G,u_1);
figure
plot(u_1,'r')
subplot(2,1,1)
plot(yu,'b')
subplot(2,1,2)
plot(yg,'r-*')
hold on
plot(phi_G*teta_G(:,i),'black--')
y_hat_G= phi_G*teta_G;
% error = sum((y-y_hat(:,length(y))).^2)





    
    
    
    
    















