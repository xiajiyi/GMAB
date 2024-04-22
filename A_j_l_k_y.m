function A_jlk= A_j_l_k_y(f11,a,b,L,c,Delta,y,m0)
A_jlk=zeros(m0,L+1,m0,length(y));
for j=0:L
    for j1=1:length(y)
        A_jlk(:,j+1,:,j1)=real(f11(:,j+1,:).*exp((-1).*1i.*(j.*pi/(b-a))*(a+c*Delta-y(j1))));
    end
end
end

