for i=1:2^6
   X(i) = (i-1)/2^6; 
end
Y=f(X);
x(1) = 0; x(1000) = 1;
for i=2:999
    x(i) = x(i-1)+1/999;
end
l = L(x,X,Y);
y=f(x);  
%误差图
for i=1:1000
    err(i)=abs(l(i)-y(i));
end
figure
semilogy(x,err);

function y=L(x,X,Y)
    n = length(X);
    m = length(x);
    for i=1:m
       y(i)=0;
       for k=1:n
          l = ((-1)^(k-1))/2^6*sin(2^6*pi*x(i))*cot(pi*(x(i)-X(k)));
          y(i) = y(i) + Y(k)*l;
       end
    end
end

function Y=f(X)
    n = length(X);
    for i=1:n
       Y(i) = sin(2*pi*X(i))*exp(cos(2*pi*X(i)));
    end
end