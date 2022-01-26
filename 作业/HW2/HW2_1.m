%%
%误差图
[x,s] = Spline(4);
n = [1:2000];
for i=1:2000
    err(i)=abs(s(i)-f(x(i)));
end
figure
semilogy(n,err);
%%
%最大误差图
for n=4:10
    [x,s] = Spline(n);
    for i=1:2000
    err(i)=abs(s(i)-f(x(i)));
    end
    err_max(n-3) = max(err);
end
n = [2^4 2^5 2^6 2^7 2^8 2^9 2^10];
figure
loglog(n,err_max);
%%
function [x,s] = Spline(n)
    k=2^n;
    h=2/k;
    X(1) = -1; X(k+1) = 1;
    for i=2:k+1
        X(i) = X(i-1)+h;
    end
    for i=1:k+1
        Y(i) = f(X(i));
    end
    x(1) = -1; x(2000) = 1;
    for i=2:1999
        x(i) = x(i-1)+2/1999;
    end
    %最后一个参数选择边界类型
    s = threesimple_1(X,Y,x,k,1);
end
%第二类边界条件
function [M,h]=threesimple_2(Y,k)
    h=2/k;
    lambda = h/(h+h);
    miu = 1-lambda;
    for i=2:k
        d(i) = (6/(h+h))*((Y(i+1)-Y(i))/h-(Y(i)-Y(i-1))/h);
    end
    d(1)=(6*(Y(2)-Y(1)))/(h*h);
    d(k+1)= (-6*(Y(k+1)-Y(k)))/(h*h);
    A = diag(repmat(2,1,k+1))+diag(repmat(miu,1,k),-1)+ ...
        diag(repmat(lambda,1,k),1);
    A(1,2)=1;   A(k+1,k)=1;
    D = d.'; 
    M=A\D;
end
%第三类边界条件
function [M,h]=threesimple_3(Y,k)
    h=2/k;
    lambda = h/(h+h);
    miu = 1-lambda;
    for i=2:k
        d(i-1)=(6/(h+h))*((Y(i+1)-Y(i))/h-(Y(i)-Y(i-1))/h);
    end
    d(k)= (6/(h+h))*((Y(2)-Y(1))/h-(Y(1)-Y(k))/h);
    A = diag(repmat(2,1,k))+diag(repmat(miu,1,k-1),-1)+ ...
        diag(repmat(lambda,1,k-1),1);
    A(1,k)=miu;   A(k,1)=lambda; A(k,k-1)=miu; A(k,k)=2;
    D = d.'; 
    N=A\D;
    M=[N(k),N.'].';
end
%插值函数
function s = threesimple_1(X,Y,x,n,flag)
    if flag == 1
        [M,h] = threesimple_2(Y,n);
    else
        [M,h] = threesimple_3(Y,n);
    end
    n=length(X); k=length(x);    
    for t=1:k
       for i=1:n-1
          if ((x(t)>=X(i)) && (x(t)<=X(i+1)))
             p1=(M(i)*(X(i+1)-x(t))^3+ ...
                 M(i+1)*(x(t)-X(i))^3)/(6*h);
             p2=( Y(i)*(X(i+1)-x(t))+ ...
                 Y(i+1)*(x(t)-X(i)) )/h;
             p3=( M(i)*(X(i+1)-x(t))+ ...
                 M(i+1)*(x(t)-X(i)) )*h/6;
             s(t)=p1+p2-p3; 
             break;
          else
             s(t)=0;
          end
       end
   end
end
%被插函数
function y = f(t)
    y = sin(4*t^2) + (sin(4*t))^2;
end