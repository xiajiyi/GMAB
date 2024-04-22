function AuxFun_2 = AuxFun_2(a,b,L,d1,d2)
%INT_1 int_d1^d2 sin(l*pi*(x-a)/(b-a)),x)
f=zeros(1,L+1);
for l=0:L
    if l==0
        f(1)=0;
    else
        f(l+1)=(b-a)/(l*pi)*(cos(l*pi*(d1-a)/(b-a))-cos(l*pi*(d2-a)/(b-a)));
    end
end
AuxFun_2 = f;
end
