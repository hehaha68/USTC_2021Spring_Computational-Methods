exact=7.9549265210128452;
N=20; %分割区间
err=zeros(1,N);
m=linspace(1,N,N);

for i=1:N
    h=2*pi/m(i);
    x=linspace(-pi,pi,m(i)+1);
    ans=0.5*f(x(1))+0.5*f(x(m(i)+1));
    for j=2:m(i)
        ans=ans+f(x(j));
    end
    ans=ans*h;
    err(i)=abs(ans-exact);
end

semilogy(m,err);

function y=f(x)
    y=exp(cos(x));
end