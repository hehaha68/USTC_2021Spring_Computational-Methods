x(1) = -1; x(2000) = 1;
for i=2:1999
    x(i) = x(i-1)+2/1999;
end
y = f(x);

for k=2:7
    %X = Xj(k);
    X = RandXj(k)
    Y = f(X);
    [A,N]=Newton(x,X,Y);
    err=zeros(1,2000);
    for i=1:2000
    err(i)=abs(N(i)-y(i));
    end
    err_max(k-1) = max(err);
end
%误差图
n = [2^2 2^3 2^4 2^5 2^6 2^7];
figure
semilogy(n,err_max);
%牛顿插值
function [A,N]=Newton(x,X,Y)
    A = Difference(X,Y);
    n = length(X);
    m = length(x);
    for i = 1:m
        N(i)=0;
        for j = 1:n
            a(i)=1;
            for k=2:j
                a(i) = a(i)*(x(i)-X(k-1));
            end
            N(i)=N(i)+a(i)*A(j,j+1);
        end
    end
end
%差商表
function A = Difference(X,Y) 
    n = length(X); 
    A = zeros(n,n+1);
    A(:,1) = X'; 
    A(:,2) = Y'; 
    for j = 3:n+1     
        for i = j-1:n         
            A(i,j)=(A(i,j-1)-A(i-1,j-1)) ...
                /(A(i,1)-A(i-j+2,1));     
        end
    end
end
%倒序插点
function x = Xj(k)
    n = 2^k;
    j = n;
    for i=1:n+1
        x(i) = cos(j*pi/n);
        j = j-1;
    end
end
%随机插点
function x = RandXj(k)
    m = 2^k;
    rng(22);
    n = randperm(m+1);
    for i=1:m+1
        x(i) = cos((n(i)-1)*pi/(m+1));
    end
end

function y = f(x)
    n = length(x); 
    for i=1:n
        y(i)=1/(1+25*x(i)*x(i));
    end
end
