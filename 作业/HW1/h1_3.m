fprintf('\n第二问\n');
A=[-148 -105 -83 -67;488 343 269 216; ...
    -382 -268 -210 -170;50 38 32 29];
[a1,a2,v1,v2] = PM(A);
show(a1,a2,v1,v2)
[a1,a2,v1,v2] = PM(-A);
show(a1,a2,v1,v2)
fprintf('\n第三问\n');
A=[222 580 584 786;-82 -211 -208 -288; ...
    37 98 101 132;-30 -82 -88 -109];
[a1,a2,v1,v2] = PM(A);
show(a1,a2,v1,v2)
fprintf('\n第四问\n');
rng(2);
A=rand(100,100);
[a1,a2,v1,v2] = IPM(A,0.8-0.6i);
show(a1,a2,v1,v2)
%为了保证精度，采用两个独立for循环分别处理只有
%一个模最大的特征值和两个互为相反数的特征值的情
%况,所以预先认为题目为可解类型
%幂法
function [a1,a2,v1,v2] = PM(A)
    q0 = ones(size(A, 1), 1);   
    tol = 1e-15;
    %一个正值或负值 
    for k = 1:1000
        q0_bar = q0/defnorm(q0);%规范化
        q1 = A*q0_bar;
        q1_bar = q1/defnorm(q1);
        %defnorm返回模最大的分量值
        if(norm(q1_bar-q0_bar,inf) < tol)|| ...
                (norm(q1_bar+q0_bar,inf)<tol)
            a1 = defnorm(q1); a2 = a1;
            v1 = q1_bar;
            v1 = v1/norm(v1,2); v2 = v1;
            return;
        end
        q0 = q1;
    end
    %互为相反数
    q0 = ones(size(A,1),1); 
    for k = 1:1000
        q0_bar = q0/defnorm(q0);
        q1 = A*q0_bar;
        q1_bar = q1/defnorm(q1);
        q2 = A*q1_bar;
        q2_bar = q2/defnorm(q2);
        if (norm(q2_bar-q0_bar,inf) < tol)
            a1 = sqrt(defnorm(A*A*q0_bar) ...
                /defnorm(q0_bar));	
            a2 = -a1;
            v1 = a1*q0_bar+q1; v1 = v1/norm(v1,2);
            v2 = a1*q0_bar-q1; v2 = v2/norm(v2,2);
            return;
        end
        disp(sqrt(defnorm(A*A*q0_bar) ...
                /defnorm(q0_bar)));fprintf('\n');
        q0 = q1;
    end
    fprintf('超过最大迭代次数\n');
end

%反幂法
function [a1,a2,v1,v2] = IPM(A,p)
    I=eye(size(A));
    B=A-p*I;
    q0 = ones(size(B, 1), 1);   
    tol = 1e-13;
    %一个正值或负值
    for k = 1:100
        q0_bar = q0/defnorm(q0);
        q1 = LU(B,q0_bar);
        q1_bar = q1/defnorm(q1);
        %defnorm返回模最大的分量值
        if(norm(q1_bar-q0_bar,inf) < tol)|| ...
                (norm(q1_bar+q0_bar,inf)<tol)
            a1 = defnorm(q1); 
            a1 = p+1/a1;%反幂法求得原特征值
            a2 = a1;
            v1 = q1_bar;
            v1 = v1/norm(v1,2); v2 = v1;
            return;
        end
        disp(p+1/defnorm(q1));fprintf('\n');
        q0 = q1;
    end
     %互为相反数
    q0 = ones(size(A, 1), 1); 
    for k = 1:100
        q0_bar = q0/defnorm(q0);
        q1 = LU(B,q0_bar);
        q1_bar = q1/defnorm(q1);
        q2 = LU(B,q1_bar);
        q2_bar = q2/defnorm(q2);
        if (norm(q2_bar - q0_bar, inf) < tol)
            a1 = sqrt(defnorm(LU(B*B,q0_bar)) ...
                /defnorm(q0_bar));	
            a1 = p+1/a1;%反幂法求得原特征值
            a2 = -a1;
            v1 = a1*q0_bar+q1; v1 = v1/norm(v1,2);
            v2 = a1*q0_bar-q1; v2 = v2/norm(v2,2);
            return;
        end
        q0 = q1;
    end
    fprintf('超过最大迭代次数\n');
end

%考虑复数的带符号的模最大分量
function x=defnorm(v)
    [~,n]=max(abs(v));
    x=v(n);
end

%LU分解法解方程
function [x] = LU(A,b) 
    [n,m] = size(A);
    for j = 1:m    
        U(1,j) = A(1,j);
    end     
    L(1,1) = 1;
    for i = 2:n    
        L(i,1) = A(i,1)/A(1,1); 
    end 
    for i = 2:n  
        for j = i:m 
            sum = 0;
            for k = 1:i-1 
                sum = sum+L(i,k)*U(k,j);
            end 
            U(i,j) = A(i,j) - sum;
            sum = 0;
            for p = 1:i-1 
                sum = sum+L(j,p)*U(p,i);
            end 
            L(j,i) = (A(j,i)-sum)/U(i,i);
        end 
    end 
    %解Ly=b
    x = zeros(n,1);
    y(1,1) = b(1,1);
    for i = 2:n 
        sum = 0;
        for j = 1:i-1
            sum = sum+L(i,j)*y(j,1);
        end 
        y(i,1) = b(i) - sum;
    end 
    %解Ux=y
    x(n,1) = y(n,1)/U(n,m);
    for i=n-1:-1:1 
        sum = 0;
        for j = i+1:n 
            sum = sum+U(i,j)*x(j,1);
        end 
        x(i,1) = (y(i,1)-sum)/U(i,i);
    end 
end

%考虑复数的格式化结果打印
function show(a1,a2,v1,v2)
    if(a1==a2)
        fprintf('The eigenvalue is:\n');
        if(imag(a1)==0)
            fprintf('%18.16f\n\n',a1);
        elseif(imag(a1)>0)
            fprintf('%18.16f+%18.16fi\n\n' ...
                ,real(a1),imag(a1));
        else
            fprintf('%18.16f%18.16fi\n\n' ...
                ,real(a1),imag(a1));
        end
        fprintf('The eigenvector is:\n ');
        v1 = regexprep(num2str(v1'),'\s*',',');
        disp(['(',v1,')^T']);
    else
        fprintf('The eigenvalues are:\n');
        if(imag(a1)==0)
            fprintf('%18.16f',a1);
            fprintf('  and  %18.16f\n\n',a2);
        elseif(imag(a1)>0)
            fprintf('%18.16f+%18.16fi',real(a1),imag(a1));
            fprintf('  and  %18.16f%18.16fi\n\n' ...
                ,real(a2),imag(a2));
        else
            fprintf('%18.16f%18.16fi',real(a1),imag(a1));
            fprintf('  and  %18.16f+%18.16fi\n\n' ...
                ,real(a2),imag(a2));
        end
        fprintf('The eigenvectors are:\n')
        v1 = regexprep(num2str(v1'),'\s*',',');
        v2 = regexprep(num2str(v2'),'\s*',',');
        disp(['(',v1,')^T']);
        disp(['(',v2,')^T']);
    end
end