clear
tic
%%%-----------参数设置-------------%%%
F0=100;G0=100;Delta=1;c=0.1;N=10;T=N*Delta;sig_x=0.1;kappa=0.1;
%%%-----------CTMC参数-------------%%%
m0=20;
%%%-----------CIR------------------%%%
% alpha1=2;sig_r=0.01;theta=0.02;r0=0.3;
alpha1=2;sig_r=0.1;theta=0.035;r0=0.03;
% alpha1=2;sig_r=0.001;theta=0.3;r0=0.3;
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
RHO=[-0.8,-0.6,-0.4,-0.2,0,0.2,0.4,0.6,0.8];
for jjjj=1:length(RHO)
    jjjj
    rho=RHO(jjjj);
L=2^8;ll=0:L;LL=10;hh=0.01;
fh=fourierX_PE_XX_CIR(sig_x,(-2:2).*1i.*hh,Q,m0,R,Delta,alpha1,theta,sig_r,rho,r0);
c11=(fh(:,2,:)-fh(:,4,:))./2./hh;
c21=(fh(:,2,:)-2.*fh(:,3,:)+fh(:,4,:))./hh.^2;
c31=(fh(:,1,:)-2.*fh(:,2,:)+2.*fh(:,4,:)-fh(:,5,:))./2./hh.^3;
c41=(fh(:,1,:)-4.*fh(:,2,:)+6.*fh(:,3,:)-4.*fh(:,4,:)+fh(:,5,:))./hh.^4;

k1=c11;
k2=c21-c11.^2;
k4=c41-4.*c31.*c11-3.*c21.^2+12.*c21.*c11.^2-6.*c11.^4;
a2=squeeze(k1.*T-LL.*(abs(k2).*T+sqrt(abs(k4).*T)).^(1/2));
b2=squeeze(k1.*T+LL.*(abs(k2).*T+sqrt(abs(k4).*T)).^(1/2));
a12=min(min(a2));
b12=max(max(b2));


Q1=300;
qq=0:Q1-1;
Delta_z=(b12-a12)/Q1;
z=a12+(qq+1/2)*Delta_z;
z=z';
y=z;

H=max(exp(z),1);
H=repmat(H,1,m0);

fp=fourierX_PE_XX_CIR(sig_x,ll.*pi/(b12-a12),Q,m0,R,Delta,alpha1,theta,sig_r,rho,r0);%m0*K*m0


A_jlk=A_j_l_k_y(fp,a12,b12,L,c,Delta,y,m0);
A_jlk=A_jlk.*sqrt(Q1/2);
A_jlk(:,1,:,:)=A_jlk(:,1,:,:)/sqrt(2);

A_jlk_0=A_j_l_k_y(fp,a12,b12,L,c,Delta,0,m0);
A_jlk_0=A_jlk_0.*sqrt(Q1/2);
A_jlk_0(:,1,:)=A_jlk_0(:,1,:)/sqrt(2);

% 
% for n=19:-1:16
%     for l=1:m0
%         X=2/(b12-a12)*H(:,l)*Delta_z;
%         Y=dct(X);
%         Y=Y(1:L);
%         Z(:,l)=Y;
%         for k=1:L
%             C(:,k,l,:)=Z(k,l)*A_jlk(:,k,l,:);
%         end
%     end
%     CC=squeeze(sum(squeeze(sum(C,2)),2));
%     CC=CC.';
%     S=exp(-kappa*(N-n)*Delta+z);
%     S=repmat(S,1,m0);
%     for j=1:m0
%         H(:,j)=max(CC(:,j),S(:,j));
%     end
% end
% 
% for l=1:m0
%     X=2/(b12-a12)*H(:,l)*Delta_z;
%     Y=dct(X);
%     Y=Y(1:L);
%     Z(:,l)=Y;
%     for k=1:L
%         C1(:,k,l)=Z(k,l)*A_jlk_0(:,k,l);
%     end
% end
% CC=squeeze(sum(squeeze(sum(C1,2)),2));
% CC=CC.';
% CC=max(1,exp(y))*CC;
% S=exp(-kappa*(N-15)*Delta+z);
% S=repmat(S,1,m0);
% for j=1:m0
%     H(:,j)=max(CC(:,j),S(:,j));
% end
% 
% for n=14:-1:11
%     for l=1:m0
%         X=2/(b12-a12)*H(:,l)*Delta_z;
%         Y=dct(X);
%         Y=Y(1:L);
%         Z(:,l)=Y;
%         for k=1:L
%             C(:,k,l,:)=Z(k,l)*A_jlk(:,k,l,:);
%         end
%     end
%     CC=squeeze(sum(squeeze(sum(C,2)),2));
%     CC=CC.';
%     S=exp(-kappa*(N-n)*Delta+z);
%     S=repmat(S,1,m0);
%     for j=1:m0
%         H(:,j)=max(CC(:,j),S(:,j));
%     end
% end
% 
% for l=1:m0
%     X=2/(b12-a12)*H(:,l)*Delta_z;
%     Y=dct(X);
%     Y=Y(1:L);
%     Z(:,l)=Y;
%     for k=1:L
%         C1(:,k,l)=Z(k,l)*A_jlk_0(:,k,l);
%     end
% end
% CC=squeeze(sum(squeeze(sum(C1,2)),2));
% CC=CC.';
% CC=max(1,exp(y))*CC;
% S=exp(-kappa*(N-10)*Delta+z);
% S=repmat(S,1,m0);
% for j=1:m0
%     H(:,j)=max(CC(:,j),S(:,j));
% end

for n=9:-1:6
    for l=1:m0
        X=2/(b12-a12)*H(:,l)*Delta_z;
        Y=dct(X);
        Y=Y(1:L);
        Z(:,l)=Y;
        for k=1:L
            C(:,k,l,:)=Z(k,l)*A_jlk(:,k,l,:);
        end
    end
    CC=squeeze(sum(squeeze(sum(C,2)),2));
    CC=CC.';
    S=exp(-kappa*(N-n)*Delta+z);
    S=repmat(S,1,m0);
    for j=1:m0
        H(:,j)=max(CC(:,j),S(:,j));
    end
end

for l=1:m0
    X=2/(b12-a12)*H(:,l)*Delta_z;
    Y=dct(X);
    Y=Y(1:L);
    Z(:,l)=Y;
    for k=1:L
        C1(:,k,l)=Z(k,l)*A_jlk_0(:,k,l);
    end
end
CC=squeeze(sum(squeeze(sum(C1,2)),2));
CC=CC.';
CC=max(1,exp(y))*CC;
S=exp(-kappa*(N-5)*Delta+z);
S=repmat(S,1,m0);
for j=1:m0
    H(:,j)=max(CC(:,j),S(:,j));
end

for n=4:-1:1
    for l=1:m0
        X=2/(b12-a12)*H(:,l)*Delta_z;
        Y=dct(X);
        Y=Y(1:L);
        Z(:,l)=Y;
        for k=1:L
            C(:,k,l,:)=Z(k,l)*A_jlk(:,k,l,:);
        end
    end
    CC=squeeze(sum(squeeze(sum(C,2)),2));
    CC=CC.';
    S=exp(-kappa*(N-n)*Delta+z);
    S=repmat(S,1,m0);
    for j=1:m0
        H(:,j)=max(CC(:,j),S(:,j));
    end
end

for l=1:m0
    X=2/(b12-a12)*H(:,l)*Delta_z;
    Y=dct(X);
    Y=Y(1:L);
    Z(:,l)=Y;
    for k=1:L
        C(:,k,l,:)=Z(k,l)*A_jlk(:,k,l,:);
    end
end
CC=squeeze(sum(squeeze(sum(C,2)),2));
CC=CC.';

j0=floor(-a12*Q1/(b12-a12)-1/2);
res(jjjj)=CC(j0,alpha0)*F0;
end
hold on
plot(RHO,res)
toc
