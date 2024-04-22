clear
tic
%%%-----------参数设置-------------%%%
F0=100;G0=100;Delta=5;c=0.1;NNN=2;T=Delta*NNN;sig_x=0.1;
mu_J=-0.1;sigma_J=0.1;lambda_J=0.2;%MJD参数
%%%-----------CTMC参数-------------%%%
m0=20;
%%%-----------CIR------------------%%%
% alpha1=2;sig_r=0.001;theta=0.3;r0=0.3;
alpha1=2;sig_r=0.1;theta=0.035;r0=0.03;
%%%-------计算网格边界方法1---------%%%
nu_bar=0.00001;gamma=1;tt=T/2;
mu_bar=r0.*exp(-alpha1.*tt)+theta.*(1-exp(-alpha1.*tt));
var_bar=sig_r^2/alpha1.*r0.*(exp(-alpha1.*tt)-exp(-2.*alpha1.*tt))+theta.*sig_r^2/(2.*alpha1).*(1-exp(-alpha1.*tt)).^2;
R=zeros(1,m0);
R(1)=max(nu_bar,mu_bar-gamma.*sqrt(var_bar));
R(m0)=mu_bar+gamma.*sqrt(var_bar);
%%%-----------计算R_j--------------%%%
alphabar=(R(m0)-R(1))/2;
c1=asinh((R(1)-r0)./alphabar);
c2=asinh((R(m0)-r0)./alphabar);
R(1,2:m0-1)=r0+alphabar.*sinh(c2.*(2:m0-1)./m0+c1.*(1-(2:m0-1)./m0));
[~,Index]=min(abs(R-r0));
R=R+r0-R(Index);
kk=R(2:m0)-R(1:m0-1);
alpha0=find(abs(R-r0)<=0.0000001);
%%%----------计算转移强度矩阵Q------%%%
mu_hat=alpha1.*(theta-R);%mu
sigma_hat=sig_r.*sqrt(R);%
Q=zeros(m0);
for j=2:m0-1
    Q(j,j-1)=((sigma_hat(j))^2-mu_hat(j)*kk(j))./(kk(j-1)*(kk(j-1)+kk(j)));
    Q(j,j+1)=((sigma_hat(j))^2+mu_hat(j)*kk(j-1))./(kk(j)*(kk(j-1)+kk(j)));
end
Q(1,2)=((sigma_hat(1))^2+mu_hat(1)*kk(1))./(kk(1)*(kk(1)+kk(1)));
Q(m0,m0-1)=((sigma_hat(m0-1))^2-mu_hat(m0-1)*kk(m0-1))./(kk(m0-1)*(kk(m0-1)+kk(m0-1)));
Q(1,1)=-Q(1,2);
Q(m0,m0)=-Q(m0,m0-1);
for j=2:m0-1
    Q(j,j)=-Q(j,j-1)-Q(j,j+1);
end
%%%-------------COS方法----------%%%
rho=0.01;
    
L=2^8;ll=0:L;LL=8;hh=0.01;
fh=MJD_fourierX_PE_XX_2_CIR(mu_J,sigma_J,lambda_J,sig_x,(-2:2).*1i.*hh,Q,m0,R,Delta,alpha1,theta,sig_r,rho,r0);
c11=(fh(:,2)-fh(:,4))./2./hh;
c21=(fh(:,2)-2.*fh(:,3)+fh(:,4))./hh.^2;
c31=(fh(:,1)-2.*fh(:,2)+2.*fh(:,4)-fh(:,5))./2./hh.^3;
c41=(fh(:,1)-4.*fh(:,2)+6.*fh(:,3)-4.*fh(:,4)+fh(:,5))./hh.^4;

k1=c11;
k2=c21-c11.^2;
k4=c41-4.*c31.*c11-3.*c21.^2+12.*c21.*c11.^2-6.*c11.^4;
a2=squeeze(k1.*T-LL.*(abs(k2).*T+sqrt(abs(k4).*T)).^(1/2));
b2=squeeze(k1.*T+LL.*(abs(k2).*T+sqrt(abs(k4).*T)).^(1/2));
a12=min(min(a2));
b12=max(max(b2));


%%%-----------COS-------------%%%
f11=MJD_fourierX_PE_XX_2_CIR(mu_J,sigma_J,lambda_J,sig_x,ll.*pi/(b12-a12),Q,m0,R,Delta,alpha1,theta,sig_r,rho,r0);%m0*K
A_jlk=A_j_l_k_1(f11,a12,b12,L);
A_jlk(:,1)=A_jlk(:,1)/2;
XS_1=-exp(-c*Delta).*AuxFun_3(a12,b12,L,a12,c*Delta)+AuxFun_1(a12,b12,L,a12,c*Delta);
E_j_1=A_jlk*XS_1';
E_j_2=exp(-c*Delta).*MJD_fourierX_PE_XX_1_CIR(mu_J,sigma_J,lambda_J,sig_x,Q,R,Delta,alpha1,theta,sig_r,rho,r0);
E_j_2=sum(E_j_2,2);
E_j=E_j_1+E_j_2;
E1=F0*(E_j_1(alpha0).^NNN)
E2=F0*(E_j_2(alpha0).^NNN)
E=E_j(alpha0);
V0=F0*E.^NNN
toc

