error=1e-15; %误差上限
x0=[1 0 0 0]';%x0是初始特征向量
%第二问
A=[-148 -105 -83 -67;
    488 343 269 216;
    -382 -268 -210 -170;
    50 38 32 29];
[x1,x2,m1,m2,tag]=eigen(A,x0,error);
show(x1,x2,m1,m2,tag);
[x1,x2,m1,m2,tag]=eigen(-A,x0,error);%测试-A的情况
show(x1,x2,m1,m2,tag);
%第三问
A=[222 580 584 786;
    -82 -211 -208 -288;
    37 98 101 132;
    -30 -82 -88 -109];
[x1,x2,m1,m2,tag]=eigen(A,x0,error);
show(x1,x2,m1,m2,tag);
% 第四问
rng(2);
A=rand(100,100);
x0=ones(100,1);%重新选择初始向量
[x1,x2,m1,m2,tag]=ipow(A,x0,error,0.8-0.6i);%反幂法
show(x1,x2,m1,m2,tag);
function [x1,x2,m1,m2,tag]=eigen(A,x0,error)%求模最大的特征值与特征向量
    m1=0; m2=0;n=0;
    x1=xbar(x0);%无穷范数规范化
    tag=1; %tag=1表示是否有两个模相等大小相反的特征值，默认为1
    for i=1:1000
        n=n+1;
        x2=A*x1;
        if norm(xbar(x2)-x1)<error 
            m1=myinf(x2);m2=m1;tag=0;
            break;
        end
        x1=xbar(x2);
    end
    if tag==0
        x1=xbar(x1);x2=xbar(x2);
        return;
    end
    x3=A*x2;  %运行到这里说明确实有模一样大的特征值
    m1=sqrt(myinf(x3)/myinf(x1));m2=-m1;
    t=x3+m1*x2;x2=x3-m1*x2;x1=t;
    % 无穷范数规范化
    x1=xbar(x1);x2=xbar(x2);
end
function a=myinf(x) %得到向量x的“无穷范数”（这里规定的无穷范数可以是复数和负数）
    [~,l]=max(abs(x));
    a=x(l);
end
function x=xbar(x) %将向量x规范化（强制要求最大的分量变为1）
    x=x/myinf(x);
end
function [x1,x2,m1,m2,tag]=ipow(A,x0,error,d) %反幂法，d是平移距离
    A = A-d*eye(size(A)); %平移变换
    [x1,x2,m1,m2,tag]=LUeigen(A,x0,error);%利用LU分解迭代
    if tag==0  %将特征值还原
         m1=d+1/m1;
         m2=m1;
    else
        m1=d+1/m1;
        m2=d+1/m2;
    end
end
function show(x1,x2,m1,m2,tag)     %将特征值和特征向量显示出来
    if tag==0
        fprintf('有一个最大的特征值\n');
        fprintf('特征值为:\n');
        printev(m1);
        fprintf('特征向量为：\n');
        disp(x1);
    else
        fprintf('有两个最大的特征值\n');
        fprintf('特征值为:\n');
        printev(m1);
        printev(m2);
        fprintf('特征向量分别为：\n');
        disp(x1);
        disp(x2);
    end
end
function printev(m)%显示特征值，包括特征值为复数的情况
    if imag(m)==0
        fprintf('%.10f\n',m);
    elseif imag(m)>0
        fprintf('%.10f+%.10fi\n',m,imag(m));
    else
        fprintf('%.10f%.10fi\n',m,imag(m));
    end
end
function [L,U]=myLU(A)%进行LU分解
    n=size(A);n=n(1); %n代表矩阵A的阶数
    L=ones(n,n);      %初始化下三角阵        
    for j=1:n-1       %LU分解定义，这里A是形参
        for i=j+1:n
            m=A(i,j)/A(j,j);
            L(i,j)=m;  
            for k=j:n
                A(i,k)=A(i,k)-m*A(j,k);
            end
        end
    end
    U=A;
end 
function [x1,x2,m1,m2,tag]=LUeigen(A,x0,error)%与函数eigen几乎没变，只是迭代部分用LU方法了
    m1=0; m2=0;n=0;
    x1=xbar(x0);%规范化
    tag=1;
    [L,U]=myLU(A);   %LU分解
    for i=1:1000
        n=n+1;
        x2=iteration(L,U,x1);
        if norm(xbar(x2)-x1)<50*error %此时只有一个模最大的特征值;多次测试发现精度只能达到1e-13
            m1=myinf(x2);
            m2=m1;
            tag=0;
            break;
        end
        x1=xbar(x2);
    end
    if tag==0
        x1=xbar(x1);
        x2=xbar(x2);
        return;
    end
    fprintf('!\n');%迭代次数超过上限，有精确度不足的风险
    
    x3=iteration(L,U,x2);  
    m1=sqrt(myinf(x3)/myinf(x1));%获得特征值
    m2=-m1;
    
    t=x3+m1*x2;
    x2=x3-m1*x2;
    x1=t;
    %无穷范数规范化
    x1=xbar(x1);
    x2=xbar(x2);
end
function x=iteration(L,U,x0)  %利用LU分解迭代，x= A_逆 * x0  
    y=x0;n=length(x0);
    y(1)=x0(1);
    for i=2:n    
        for j=1:i-1
            x0(i)=x0(i)-L(i,j)*y(j);
        end
        y(i)=x0(i);
    end

    x(n)=y(n)/U(n,n);
    for i=(n-1):-1:1
        for j=n:-1:i+1
            y(i)=y(i)-U(i,j)*x(j);
        end
        x(i)=y(i)/U(i,i);
    end
    x=x-2i.*imag(x);
    x=x';
end 