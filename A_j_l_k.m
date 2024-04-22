function A_jlk= A_j_l_k(f11,a,b,L)
A_jlk=zeros(size(f11));
for j=0:L
    A_jlk(:,j+1,:)=2/(b-a).*real(f11(:,j+1,:).*exp((-1).*1i.*(j.*pi/(b-a))*a));
end
end

