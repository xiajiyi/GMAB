function AuxFun_3 = AuxFun_3(a,b,L,d1,d2)
%INT_1 int_d1^d2 e^x cos(l*pi*(x-a)/(b-a)),x)
f=zeros(1,L+1);
for l=0:L
    f(l+1)=1/(1+(l*pi/(b-a))^2)*((exp(d2)*cos(l*pi*(d2-a)/(b-a))-exp(d1)*cos(l*pi*(d1-a)/(b-a)))...
        +l*pi/(b-a)*(exp(d2)*sin(l*pi*(d2-a)/(b-a))-exp(d1)*sin(l*pi*(d1-a)/(b-a))));
end
AuxFun_3 = f;
end
