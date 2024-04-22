function PSI=fourierX_PE_XX_3_CIR(sig_x,w,Q,m0,R,Delta,alpha1,lam_R,sig_R,rho,r0)
f2_R=alpha1.*(lam_R-R).*sig_x./(sig_R.*sqrt(R))-1./4.*sig_R.*sig_x./(sqrt(R));
f3_R=2.*rho.*sig_x./sig_R.*(sqrt(R)-sqrt(r0));
phi_f1_fis=(1i*w.'-1).*R-rho*1i.*w.'.*f2_R;%w*m0
j1=length(w);
KKK=zeros(m0,j1,m0);
for j=1:j1
    KKK(:,j,:)=expm(Delta.*(Q+diag(phi_f1_fis(j,:))));
end
phi_X=exp(-1i.*w.'*Delta*sig_x^2/2-sig_x^2*(1-rho^2).*Delta.*w.'.^2/2+1i.*w.'*f3_R);
PSI=zeros(m0,j1,m0);
for j=1:m0
    PSI(j,:,:)=phi_X.*squeeze(KKK(j,:,:));
end
end