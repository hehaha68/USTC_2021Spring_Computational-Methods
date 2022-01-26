A = [2 -1 0 0 0 0 0 0 0 0;
    -1 2 -1 0 0 0 0 0 0 0;
    0 -1 2 -1 0 0 0 0 0 0;
    0 0 -1 2 -1 0 0 0 0 0;
    0 0 0 -1 2 -1 0 0 0 0;
    0 0 0 0 -1 2 -1 0 0 0;
    0 0 0 0 0 -1 2 -1 0 0;
    0 0 0 0 0 0 -1 2 -1 0;
    0 0 0 0 0 0 0 -1 2 -1;
    0 0 0 0 0 0 0 0 -1 2 ];
b=[2 -2 2 -1 0 0 1 -2 2 -2]';
superior=600;%迭代次数上限
exact=[1 0 1 0 0 0 0 -1 0 -1]';%精确解
error=1e-15;%误差
xold=[1 1 1 1 1 1 1 1 1 1]';
tag=1;%是否优化:0代表不优化，1代表优化
[y1,x]=Jacob(A,b,xold,error,exact,superior,tag);
fprintf('jacob迭代结果为：\n');
disp(x);
[y2,x]=Gauss_Seidel(A,b,xold,error,exact,superior,tag);
fprintf('Gauss_Seidel迭代结果：\n');
disp(x);
[y3,x]=SOR(A,b,1.666,xold,error,exact,superior,tag);
fprintf('松弛因子为1.666的SOR迭代结果:\n');
disp(x);
[y4,x]=SOR(A,b,1.5,xold,error,exact,superior,tag);
fprintf('松弛因子为1.5的SOR迭代结果:\n');
disp(x);
[y5,x]=SOR(A,b,1.75,xold,error,exact,superior,tag);
fprintf('松弛因子为1.75的SOR迭代结果:\n');
disp(x);
C=compare(A,b,xold,error,exact,superior,20);
fprintf('优化前与优化后19次累计用时（秒）：\n');
disp(C);

%绘图
n1=1:length(y1);
n2=1:length(y2);
n3=1:length(y3);
n4=1:length(y4);
n5=1:length(y5);
semilogy(n1,y1);
text(400,1e-7,'Jacob');
hold on
semilogy(n2,y2);
text(203,1e-7,'Gauss Seidel');
hold on;
semilogy(n3,y3);
text(78,2e-14,'SOR 1.666');
hold on;
semilogy(n4,y4);
text(87,7e-9,'SOR 1.5');
semilogy(n5,y5);
text(85,4e-11,'SOR 1.7');
title '迭代次数与误差的关系'；
legend('Jacob','Gauss Seidel','SOR 1.666','SOR 1.5','SOR 1.7');

%Jacob方法:
function [y,xnew]=Jacob(A,b,xold,error,exact,superior,tag)
    xnew=[0 0 0 0 0 0 0 0 0 0]';
    B1= -0.5 .* eye(size(A)) *(A- 2 .* eye(size(A)));%迭代矩阵
    d1= 0.5 .* eye(size(A))*b;
    n=1;
    y=[];
    while 1
        if n<=superior
            y=[y,norm(xold-exact)];
            if tag==0  %判断是否优化
                xnew=B1*xold+d1;
            else
                for i=1:10
                    if i==1
                        xnew(1)=B1(1,2)*xold(2)+d1(i);
                    elseif i==10
                        xnew(10)=B1(10,9)*xold(9)+d1(i);
                    else
                        xnew(i)=B1(i,i-1)*xold(i-1)+B1(i,i+1)*xold(i+1)+d1(i);
                    end
                end
            end
            if norm(xnew-exact)<error
                break;
            end
            xold=xnew;
            n=n+1;
        else
            fprintf('迭代次数达到上限\n');
            break;
        end
    end

end

%Gauss-Seidel方法
function [y,xnew]=Gauss_Seidel(A,b,xold,error,exact,superior,tag)
    n=0;
    y=[];
    xnew=[0 0 0 0 0 0 0 0 0 0]';
    while 1
        if tag==0
            xnew(1)=(b(1)-A(1,2:10)*xold(2:10))/A(1,1);
            for i=2:10
                if i==10
                    xnew(10)=(b(10)-A(10,1:9)*xnew(1:9))/A(10,10);
                else
                    xnew(i)=(b(i)-A(i,1:i-1) * xnew(1:i-1)-A(i,i+1:10) * xold(i+1:10))/A(i,i);
                end
            end
        else %优化模式
            xnew(1)=(b(1)-A(1,2)*xold(2))/A(1,1);
            for i=2:10
                if i==10
                    xnew(10)=(b(10)-A(10,9)*xnew(9))/A(10,10);
                else
                    xnew(i)=(b(i)-A(i,i-1) * xnew(i-1)-A(i,i+1) * xold(i+1))/A(i,i);
                end
            end
        end
        n=n+1;
        y=[y,norm(xnew-exact)];
        if norm(xnew-xold)<error
            break;
        end
        xold=xnew;%更新解向量
        if n>superior
            fprintf('迭代次数达到上限\n');
            break;
        end
    end 

end
% SOR方法
function [y,xnew]=SOR(A,b,w,xold,error,exact,superior,tag)
    n=0;
    y=[];
    xnew=[0 0 0 0 0 0 0 0 0 0]';
    while 1
        if tag==0
            xnew(1)=(b(1)-A(1,2:10) * xold(2:10))/A(1,1);
            for i=2:10
                xnew(i)=(1-w)*xold(i) + w*(b(i)-A(i,1:i-1)*xnew(1:i-1)-A(i,i+1:10)*xold(i+1:10))/A(i,i);
            end
       
        else %优化模式
            xnew(1)=(b(1)-A(1,2) * xold(2))/A(1,1);
            for i=2:9
                xnew(i)=(1-w)*xold(i) + w*(b(i)-A(i,i-1)*xnew(i-1)-A(i,i+1)*xold(i+1))/A(i,i);
            end
            xnew(10)=(b(10)-A(10,9)*xold(9))/A(10,10);
        end
        if norm(xnew-exact)<error
            break;
        end
        n=n+1;
        y=[y,norm(xnew-exact)];
        xold=xnew;
        if n>superior
            fprintf('迭代次数达到上限\n');
            break;
        end
    end

end
function diff=compare(A,b,xold,error,exact,superior,N)
    diff=zeros(3,2);
    J0=[];GS0=[];SOR0=[];
    J1=[];GS1=[];SOR1=[];
    tag=0;%关闭优化
    for i=1:N
        tic;
        Jacob(A,b,xold,error,exact,superior,tag);
        toc;
        J0=[J0,toc];
        tic;
        Gauss_Seidel(A,b,xold,error,exact,superior,tag);
        toc;
        GS0=[GS0,toc];
        tic;
        SOR(A,b,1.65,xold,error,exact,superior,tag);
        toc;
        SOR0=[SOR0,toc];
    end
    tag=1;%启用优化
    for i=1:N
        tic;
        Jacob(A,b,xold,error,exact,superior,tag);
        toc;
        J1=[J1,toc];
        tic;
        Gauss_Seidel(A,b,xold,error,exact,superior,tag);
        toc;
        GS1=[GS1,toc];
        tic;
        SOR(A,b,1.65,xold,error,exact,superior,tag);
        toc;
        SOR1=[SOR1,toc];
    end
    for i=2:N
        diff(1,1)=diff(1,1)+J0(i);
        diff(1,2)=diff(1,2)+J1(i);
        diff(2,1)=diff(2,1)+GS0(i);
        diff(2,2)=diff(2,2)+GS1(i);
        diff(3,1)=diff(3,1)+SOR0(i);
        diff(3,2)=diff(3,2)+SOR1(i);
    end
end

