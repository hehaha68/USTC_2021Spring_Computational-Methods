syms x;
f(x)=x^3-3*x^2+2;
df(x)=diff(f(x));
fprintf('\nThe left one:\n');
NewtonIteration(f,df,-1.5);
fprintf('\nThe middle one:\n');
NewtonIteration(f,df,0.5);
fprintf('\nThe right one:\n');
NewtonIteration(f,df,3);
%牛顿迭代法
function NewtonIteration(fun,dfun,x0)
    f = fun;
    df = dfun;
    fprintf ...
        ('iter.        x             log(|err|)     order\n')
    dashString = repmat('-', 1, 50);
    fprintf('%s\n', dashString);
    x1 = x0 - f(x0)/df(x0);
    d1 = norm(x1-x0);
    k = 1;
    tol = 100;%迭代次数上限
    while d1 > eps
        if k < tol
            if k > 2 %第三次迭代开始估计收敛阶数
                fprintf ...
                    ('%d    %.15f    %.8f   %.3f\n' ...
                    ,k,x1,log(d1),log(d1)/log(d2));
            else
                fprintf('%d    %.15f    %.8f\n',k,x1,log(d1));
            end
            x0 = x1;
            x1 = x0 - f(x0)/df(x0);%迭代序列
            d2 = d1;
            d1 = norm(x1-x0);%两次迭代相差
            k = k + 1;
        else
            break;
        end
    end
    fprintf('%s\n', dashString);
    fprintf('The solutiion is: %.15f\n',x1);
end
