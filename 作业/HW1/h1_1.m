clear,clc;
A=diag(repmat(2,1,10))+diag(repmat(-1,1,9),1) ...
    +diag(repmat(-1,1,9),-1);
b=[2 -2 2 -1 0 0 1 -2 2 -2]';
x_exact=[1 0 1 0 0 0 0 -1 0 -1]';
x0=[ones(1,size(b,1))]';
%%
%误差图
[a1,b1,~]=Jacobi(A,b,x0,x_exact);
[a2,b2,~]=G_S(A,b,x0,x_exact);
[a3,b3,~]=SOR(A,b,1.5,x0,x_exact);
[a4,b4,~]=SOR(A,b,1.25,x0,x_exact);
[a5,b5,~]=SOR(A,b,1.05,x0,x_exact);
[a6,b6,~]=SOR(A,b,0.75,x0,x_exact);
[a7,b7,~]=SOR(A,b,0.55,x0,x_exact);
[a8,b8,~]=SOR(A,b,0.45,x0,x_exact);
[a9,b9,~]=SOR(A,b,1.75,x0,x_exact);
[a10,b10,~]=SOR(A,b,1.65,x0,x_exact);

n1=1:1:a1;
n2=1:1:a2;
n3=1:1:a3;
n4=1:1:a4;
n5=1:1:a5;
n6=1:1:a6;
n7=1:1:a7;
n8=1:1:a8;
n9=1:1:a9;
n10=1:1:a10;

figure
semilogy(n1,b1,n2,b2,n3,b3,n4,b4,n5,b5, ...
    n6,b6,n7,b7,n8,b8,n9,b9,n10,b10);
title('The semilogy figure');
xlabel('No. of iterations');
ylabel('Error');
text(440,8e-8,'Jacobi');
text(290,6e-10,'G-S');
text(70,3e-12,'\omega=1.5');
text(180,1e-12,'\omega=1.25')
text(300,2e-13,'\omega=1.05')
text(450,3e-10,'\omega=0.75')
text(470,4e-6,'\omega=0.55')
text(480,7e-5,'\omega=0.45')
text(130,5e-15,'\omega=1.75')
text(60,1e-15,'\omega=1.65')
%%
%稀疏矩阵改进与用时比较
for i=1:10
    [~,~,t1]=Jacobi(A,b,x0,x_exact);
    t(i)=t1;
end
T(1) = sum(t(2:10))/9;
for i=1:10
    t2=Jacobi_new(A,b,x0,x_exact);
    t(i)=t2;
end
T(2) = sum(t(2:10))/9;
for i=1:10
    [~,~,t1]=G_S(A,b,x0,x_exact);
    t(i)=t1;
end
T(3) = sum(t(2:10))/9;
for i=1:10
    t2=G_S_new(A,b,x0,x_exact);
    t(i)=t2;
end
T(4) = sum(t(2:10))/9;
for i=1:10
    [~,~,t1]=SOR(A,b,1.65,x0,x_exact);
    t(i)=t1;
end
T(5) = sum(t(2:10))/9;
for i=1:10
    t2=SOR_new(A,b,1.65,x0,x_exact);
    t(i)=t2;
end
T(6) = sum(t(2:10))/9;
for i=1:10
    [~,~,t1]=SOR(A,b,1.05,x0,x_exact);
    t(i)=t1;
end
T(7) = sum(t(2:10))/9;
for i=1:10
    t2=SOR_new(A,b,1.05,x0,x_exact);
    t(i)=t2;
end
T(8) = sum(t(2:10))/9;
for i=1:10
    [~,~,t1]=SOR(A,b,1.75,x0,x_exact);
    t(i)=t1;
end
T(9) = sum(t(2:10))/9;
for i=1:10
    t2=SOR_new(A,b,1.75,x0,x_exact);
    t(i)=t2;
end
T(10) = sum(t(2:10))/9;
for i=1:10
    [~,~,t1]=SOR(A,b,1.25,x0,x_exact);
    t(i)=t1;
end
T(11) = sum(t(2:10))/9;
for i=1:10
    t2=SOR_new(A,b,1.25,x0,x_exact);
    t(i)=t2;
end
T(12) = sum(t(2:10))/9;
for i=1:10
    [~,~,t1]=SOR(A,b,0.75,x0,x_exact);
    t(i)=t1;
end
T(13) = sum(t(2:10))/9;
for i=1:10
    t2=SOR_new(A,b,0.75,x0,x_exact);
    t(i)=t2;
end
T(14) = sum(t(2:10))/9;
for i=1:10
    [~,~,t1]=SOR(A,b,0.45,x0,x_exact);
    t(i)=t1;
end
T(15) = sum(t(2:10))/9;
for i=1:10
    t2=SOR_new(A,b,0.45,x0,x_exact);
    t(i)=t2;
end
T(16) = sum(t(2:10))/9;

dashString = repmat('-', 1, 35);
fprintf('              Oldtime     Newtime\n');
fprintf('%s\n', dashString);
fprintf('Jacobi:');
fprintf('      %fs   %fs\n',T(1),T(2));
fprintf('G-S:');
fprintf('         %fs   %fs\n',T(3),T(4));
fprintf('SOR w=1.65:');
fprintf('  %fs   %fs\n',T(5),T(6));
fprintf('SOR w=1.05:');
fprintf('  %fs   %fs\n',T(7),T(8));
fprintf('SOR w=1.75:');
fprintf('  %fs   %fs\n',T(9),T(10));
fprintf('SOR w=1.25:');
fprintf('  %fs   %fs\n',T(11),T(12));
fprintf('SOR w=0.75:');
fprintf('  %fs   %fs\n',T(13),T(14));
fprintf('SOR w=0.45:');
fprintf('  %fs   %fs\n',T(15),T(16));
fprintf('%s\n', dashString);
%%
%Jacobi
function [k,y,t]=Jacobi(A,b,x0,x_exact)
    tic;
    D = diag(diag(A));
    x_old = x0;
    x_new = [zeros(1,size(b,1))]';
    y = [norm(x_old-x_exact)];%初始误差
    k = 1;
    m = 500;%次数上限
    while norm(x_new-x_exact) > eps
        if k <= m
            x_new = D\((D-A)*x_old+b);%迭代
            y = [y,norm(x_new-x_exact)];
            k = k+1;
            x_old = x_new;
        else
            break;
        end
    end
    t = toc;
%     fprintf('The result of the Jacobi method is:\n');
%     showvector(x_new);    
end

%G-S
function [k,y,t]=G_S(A,b,x0,x_exact)
    tic;
    L = tril(A,-1);
    U = triu(A,1);
    D = diag(diag(A));
    x_old = x0;
    x_new = [zeros(1,size(b,1))]';
    y = [norm(x_old-x_exact)];%初始误差
    k = 1;
    m=500; %次数上限
    while norm(x_new-x_exact) > eps
        if k <= m         
            x_new = (D+L)\(b-U*x_old);%迭代
            y = [y,norm(x_new-x_exact)];
            x_old = x_new;
            k = k+1;
        else
            break;
        end
    end
    t = toc;
%     fprintf('The result of the Gauss-Seidel method is:\n');
%     showvector(x_new);
end

%SOR
function [k,y,t]=SOR(A,b,w,x0,x_exact)
    tic;
    L = tril(A,-1);
    U = triu(A,1);
    D = diag(diag(A));
    I=eye(size(A));
    x_old = x0;
    x_new = [zeros(1,size(b,1))]';
    y = [norm(x_old-x_exact)];%初始误差
    k = 1;
    m = 500;%次数上限
    while norm(x_new-x_exact) > eps
        if k <= m           
            x_new=(I+D\L*w)\ ...
                ((I-w*(D\U+I))*x_old+D\b*w);%迭代
            y = [y,norm(x_new-x_exact)];
            x_old = x_new;
            k = k+1;
        else
            break;
        end
    end
    t = toc;
%     fprintf ...
%     ('The result of the SOR method with w=%.2f is：\n',w);
%     showvector(x_new);
end

%Jacobi_new
function t=Jacobi_new(A,b,x0,x_exact)
    tic;
    D = diag(diag(A));
    x_old = x0;
    n = size(b,1);
    x_new = [zeros(1,n)]';
    k = 1;
    m = 500;%次数上限
    while norm(x_new-x_exact) > eps
        if k <= m
            x_new(1)=(b(1)-A(1,2)*x_old(2))/D(1,1);
            for i = 2:n-1
                x_new(i)=(b(i)-A(i,i-1)*x_old(i-1) ...
                    -A(i,i+1)*x_old(i+1))/D(i,i);
            end
            x_new(n)=(b(n)-A(n,n-1)*x_old(n-1))/D(n,n);
            k = k+1;
            x_old = x_new;
        else
            break;
        end
    end
    t = toc;
    %fprintf('The result of the Jacobi method is:\n');
    %showvector(x_new);    
end

%G-S_new
function t=G_S_new(A,b,x0,x_exact)
    tic;
    L = tril(A,-1);
    U = triu(A,1);
    D = diag(diag(A));
    n = size(b,1);
    x_old = x0;
    x_new = [zeros(1,n)]';
    k = 1;
    m = 500; %次数上限
    while norm(x_new-x_exact) > eps
        if k <= m         
            x_new(1)=(b(1)-A(1,2)*x_old(2))/D(1,1);
            for i=2:n-1
                x_new(i)=(b(i)-A(i,i-1)*x_new(i-1) ...
                    -A(i,i+1)*x_old(i+1))/D(i,i);
            end
            x_new(n)=(b(n)-A(n,n-1)*x_new(n-1))/D(n,n);
            x_old = x_new;
            k = k+1;
        else
            break;
        end
    end
    t = toc;
    %fprintf('The result of the Gauss-Seidel method is:\n');
    %showvector(x_new);
end

%SOR_new
function t=SOR_new(A,b,w,x0,x_exact)
    tic;
    L = tril(A,-1);
    U = triu(A,1);
    D = diag(diag(A));
    I = eye(size(A));
    x_old=x0;
    x_new=[zeros(1,size(b,1))]';
    n = size(b,1);
    k = 1;
    m = 500;
    while norm(x_new-x_exact) > eps
        if k <= m           
            x_new(1)=w*(b(1)-A(1,2)*x_old(2)) ...
                /D(1,1)+(1-w)*x_old(1);
            for i = 2:n-1
                x_bar=(b(i)-A(i,i-1)*x_new(i-1) ...
                    -A(i,i+1)*x_old(i+1))/D(i,i);
                x_new(i)=w*x_bar+(1-w)*x_old(i);
            end
            x_new(n)=w*(b(n)-A(n,n-1)*x_new(n-1)) ...
                /D(n,n)+(1-w)*x_old(n);
            x_old = x_new;
            k = k+1;
        else
            break;
        end
    end
    t = toc;
   % fprintf ...
   %     ('The result of the SOR method with w=%.2f is：\n',w);
   % showvector(x_new);
end

function showvector(v)%打印向量
    [n,~] = size(v);
    fprintf('(');
    fprintf('%.3f',v(1));
    for k = 2:n
        fprintf(',%.3f',v(k));
    end   
    fprintf(')^T\n\n');
end
