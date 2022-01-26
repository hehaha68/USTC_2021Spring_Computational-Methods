error=1e-15;%误差上限
%三个初始值
x0=-2;
x1=1.2;
x2=3;

xl=newton(x0,error);
fprintf('The solution of xl is %.14f\n',xl);
xm=newton(x1,error);
fprintf('The solution of xm is %.14f\n',xm);
xr=newton(x2,error);
fprintf('The solution of xr is %.14f\n',xr);

function f=f(x)
    f=x^3-3*x^2+2;
end
%f1(x)是f(x)的导数
function f=f1(x)
    f=3*x^2-6*x;
end
function resolve=newton(x0,error)
    N=100;n=0;resolve=inf;
    while 1
        %pause;
        n=n+1;
        x1=x0-f(x0)/f1(x0);
        if mod(n,2)==1
            d0=abs(x1-x0);
        else
            d1=abs(x1-x0);
        end
        fprintf('Present value: %.14f   ',x1); 
        if n<2
            fprintf('\n');
        else
            if mod(n,2)==0
                order=log10(d1)/log10(d0);
            else
                order=log10(d0)/log10(d1);
            end
            fprintf('Order:%.6f\n',order);
        end
        if abs(x1-x0)<error
            resolve=x1; 
            break;
        end
        x0=x1;%迭代
        if n>N %超出迭代次数上限，发生异常
            fprintf('ERROR\n');
            break;
        end
    end
end