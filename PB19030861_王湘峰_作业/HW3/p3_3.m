N=1000;
h=2/(N-1);
x=linspace(0,2,N);
y=zeros(1,N);
err=zeros(1,N);
x(1)=0;y(1)=0;
k1=f(x(1),y(1));
k2=f(x(1)+h/2,y(1)+h/2*k1);
y(2)=y(1)+h*k2;
err(2)=abs(Y(x(2))-y(2));
f1=f(x(1),y(1));
f2=f(x(2),y(2));
for i=3:N 
    y(i)=3/(3+4*h)*(y(i-2)+h/3*(x(i)*exp(-4*x(i))+4*f2+f1));
    f1=f2;
    f2=f(x(i),y(i));
end
plot(x,y)


function s=f(x,y)
    s=x*exp(-4*x)-4*y;
end

function y=Y(x)
    y=0.5*x^2*exp(-4*x);
end