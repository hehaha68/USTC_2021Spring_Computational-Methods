n=2^6;
x=zeros(n,1);
for i=1:n
    x(i)=(i-1)/n;
end
y=zeros(n,1);
l=zeros(n,1);
for i=1:n
    y(i)=f(x(i));
end
X=linspace(0,1,1000);
T=zeros(1000,1);
Y=T;
err=zeros(1000,1);
for i=1:1000
    T(i)=f(X(i));
    Y(i)=lagrange(X(i),x,y,n);
    err(i)=abs(T(i)-Y(i));
end
% plot(X,T);
% hold on;
% plot(X,Y);
semilogy(X,err);

function ans=lagrange(t,x,y,n)
    ans=0;
    for k=1:n
        l(k)=(-1)^(k-1)*sin(n*pi*t)*cot(pi*(t-x(k)))/n;
        ans=ans+y(k)*l(k);
    end
end
function y=f(x)
    y=sin(2*pi*x)*exp(cos(2*pi*x));
end
