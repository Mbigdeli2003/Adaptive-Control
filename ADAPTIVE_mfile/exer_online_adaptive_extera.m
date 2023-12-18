%% Exer_2_adaptive_online
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
n=100;
ind_1=0;
yg=zeros(1,n);
u_1=rand(1,n);
Po=1e4*eye(6);
teta_G=zeros(6,n);
K=zeros(6,n);
err_G=zeros(1,n)+1;
yg_test=zeros(1,n);
yy=zeros(1,n);


% P0_1=zeros
for i=4:1:max(size(u_1))
 yg(i)=-(Ag(1,2)*yg(i-1))-Ag(1,3)*yg(i-2)-(Ag(1,4)*yg(i-3))+(Bg(1,2)*u_1(i-1))+(Bg(1,3)*u_1(i-2))+(Bg(1,4)*u_1(i-3));   
 
end
for i=4:n
    phi=[yg(i-1);yg(i-2);yg(i-3);u_1(i-1);u_1(i-2);u_1(i-3)];%6x1
    K(:,i)=Po*phi*inv(1+phi'*Po*phi);%6x1
    err_G(i)=yg(i)-phi'*teta_G(:,i-1);
    yy(i)=phi'*teta_G(:,i-1);
    teta_G(:,i)=teta_G(:,i-1)+K(:,i)*err_G(i);%6x1
    Po=[eye(6)-K(:,i)*phi']*Po;%6x6    
    yg_test(i)=teta_G(1,i)*yg(i-1)+teta_G(2,i)*yg(i-2)+teta_G(3,i)*yg(i-3)+teta_G(4,i)*u_1(i-1)+teta_G(5,i)*u_1(i-2)+teta_G(6,i)*u_1(i-3);
    yg_error=abs(yg(i)-yg_test(i));
    ind_1=i;
end
for i=ind_1:n
    yy(i)=teta_G(1,ind_1)*yy(i-1)+teta_G(2,ind_1)*yy(i-2)+teta_G(3,ind_1)*yy(i-3)+teta_G(4,ind_1)*u_1(i-1)+teta_G(5,ind_1)*u_1(i-2)+teta_G(6,ind_1)*u_1(i-3);
end

figure(2)
plot(1:n,yg,'b-*')
hold on
plot(1:n,yy,'r-')






















% u_1(i-1)=inv(phi_G_1(i-1)'*phi_G_1(i-1));





    
    
    
    
    















