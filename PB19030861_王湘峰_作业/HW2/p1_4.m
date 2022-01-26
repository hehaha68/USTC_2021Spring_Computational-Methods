n=16;
x=linspace(-1,1,n+1);
y=zeros(1,n);
for i=1:n+1
    y(i)=f(x(i));
end
%M=getM(n,x,y);


l=linspace(-1,1,2000);
st=zeros(1,2000);
p=zeros(1,2000);
for i=1:2000
    st(i)=f(l(i));
end
for i=1:2000
    err(i)=abs(st(i)-myspine(l(i),n));
end
N=1:2000;
semilogy(N,err);

E=zeros(1,7);
for i=1:7
    E(i)=2^(i+3);
end
merr=zeros(1,7);
for j=1:7
    max=-1;
    for i=1:2000
        t=abs(st(i)-myspine(l(i),E(j)));
        if t>max
            max=t;
        end
    end
    merr(j)=max;
end
loglog(E,merr);

function M1=getM(n,x,y)
    %x是插值点,y是对应的函数值(数组)
    for k=1:n %计算h(i)
        h(k)=x(k+1)-x(k);
    end
    lambda=zeros(1,n);
    b=zeros(1,n);
    A=zeros(n,n);
    miu=lambda;
    for k=1:n-1 %计算miu和lambda
        lambda(k)=h(k+1)/(h(k)+h(k+1));
        miu(k)=1-lambda(k);
    end
    lambda(n)=h(1)/(h(n)+h(1));
    miu(n)=1-lambda(n);
    %disp(lambda);
    for k=1:n-1
        b(k)=6/(h(k)+h(k+1))*((y(k+2)-y(k+1))/h(k+1)-(y(k+1)-y(k))/h(k));
    end
    %b(n)=3/h(1)^2*(y(2)-2*y(n+1)-y(n));
    b(n)=6/(h(1)+h(n))*((y(2)-y(1))/h(1)-(y(1)-y(n))/h(1));
    b=b';
%   fprintf('b=');disp(b);
    A(1,n)=miu(1);A(1,1)=2;A(1,2)=lambda(1);
    A(n,1)=lambda(n);A(n,n)=2;A(n,n-1)=miu(n);
    for i=2:n-1
        A(i,i)=2;
        A(i,i-1)=miu(i);
        A(i,i+1)=lambda(i);
    end
    %disp(A);
    M=A\b;
    M1=[M(n),M']';
%   disp(M1);
end
function y=f(x)
    y=sin(4*x^2)+sin(4*x)^2;
end
function s=myspine(t,n)
    s=0;N=n;
    if t>1||t<-1
        fprintf('wrong input!\n');
        return;
    end
    x=linspace(-1,1,n+1);
    y=zeros(1,n);
    h=y;
    for i=1:n+1
        y(i)=f(x(i));
    end
    for k=1:n %计算h(i)
        h(k)=x(k+1)-x(k);
    end
    M=getM(N,x,y);
    %计算S(x)
    for i=1:n
%         if t==x(i)
%             break;
%         end
        if t>=x(i)&&t<=x(i+1)
            break;
        end
    end
        s=((x(i+1)-t)^3*M(i)+(t-x(i))^3*M(i+1))/(6*h(i))+(y(i)*(x(i+1)-t)+y(i+1)*(t-x(i)))/h(i)-h(i)/6*(M(i)*(x(i+1)-t)+M(i+1)*(t-x(i)));
    end