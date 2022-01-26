E=[4 8 16 32 64 128];
maxerr=zeros(1,6);
X=linspace(-1,1,2000);
err=zeros(1,2000);
Y=zeros(1,2000);
for i=1:2000
    Y(i)=f(X(i));
end
for j=1:6
    n=E(j);
    x=zeros(n+1,1);
    for i=0:n
        x(i+1)=cos(pi*i/n);
    end
    for i=1:2000
        err(i)=abs(Y(i)-N(X(i),x,E(j)));
    end
    maxerr(j)=max(err);
end
semilogy(E,maxerr);


function y=f(x)
    y=1/(1+25*x^2);
end
function A=getA(n,x)
    A=zeros(n+1,1);
    for i=1:n+1
        A(i)=f(x(i));
    end
    for k=2:n+1
        j=n+1;
        while j>=k
            A(j)=(A(j)-A(j-1))/(x(j)-x(j+1-k));
            j=j-1;
        end
    end
end
function y=N(t,x,n)
    A=getA(n,x);
    ans=A(1);
    h=zeros(1,n);
    for i=1:n
        h(i)=t-x(i);
    end
    for i=1:n
        prod=1;
        for j=1:i
            prod=prod*h(j);
        end
        ans=ans+prod*A(i+1);
    end
    y=ans;
end
