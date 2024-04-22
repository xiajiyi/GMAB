function PSI_1=MJD_fourierX_PE_XX_1_CIR(mu_J,sigma_J,lambda_J,sig_x,Q,R,Delta,alpha1,theta,sig_r,rho,r0)
f2_R=alpha1.*(theta-R).*sig_x./(sig_r.*sqrt(R))-1./4.*sig_r.*sig_x./(sqrt(R));
f3_R=2.*rho.*sig_x./sig_r.*(sqrt(R)-sqrt(r0));
phi_f1_fis=-rho.*f2_R;%1*m0
KKK=expm(Delta.*(Q+diag(phi_f1_fis)));
phi_X=exp(-Delta*(lambda_J*(exp(mu_J+sigma_J^2/2)-1)+sig_x^2/2)+sig_x^2*(1-rho^2).*Delta/2+f3_R+Delta*lambda_J*(exp(mu_J+sigma_J^2/2)-1));
PSI_1=phi_X.*KKK;
end
