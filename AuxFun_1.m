function AuxFun_1 = AuxFun_1(a,b,L,d1,d2)
%INT_1 int_d1^d2 cos(l*pi*(x-a)/(b-a)),x)
f=zeros(1,L+1);
for l=0:L
    if l==0
        f(1)=d2-d1;
    else
        f(l+1)=(b-a)/(l*pi)*(sin(l*pi*(d2-a)/(b-a))-sin(l*pi*(d1-a)/(b-a)));
    end
end
AuxFun_1 = f;
end
