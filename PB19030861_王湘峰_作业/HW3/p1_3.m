N=40;   %最大迭代次数
exact=cos(1.2);
err=zeros(1,N);
A=zeros(N,N);
X=zeros(1,N);
h=67;
for i=1:N
    A(i,1)=forward(h/2^(i-1));
end
for k=2:N
    for i=k:N
        A(i,k)=A(i,k-1)+(A(i,k-1)-A(i-1,k-1))/(2^(k-1)-1);
    end
end

for i=1:N
    err(i)=abs(A(i,i)-exact);
    X(i)=i;
end
% semilogy(X,err);
fprintf('%.16f\n',A(14,14));
fprintf('%.16f\n',exact);
function y=f(x)
    y=sin(x);
end
function y=forward(h)
    y=(f(1.2+h)-f(1.2))/h;
end