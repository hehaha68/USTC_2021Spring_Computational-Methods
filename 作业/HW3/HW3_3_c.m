clear;
%%
y = result(1000);
plot(y);
%%
err = zeros(1,991);
real = 1/2*exp(-4*2)*2^2;
for i=10:1000
    y = result(i);
    err(i-9) = abs(y(i+1)-real);
end
h = zeros(1,991);
for i=10:1000
    h(i-9) = 2/i;
end
figure
loglog(h,err);

%%
function y = result(n)
    h=2/n;
    x(1)=0;
    for i=2:n+1
       x(i) = x(i-1)+h; 
    end
    y = zeros(1,n+1);
    y(1)=0;
    y(2)=RK2(y(1),h,x(1));
    
    for i=3:n+1
        y(i) = MS(x,y,h,i);
    end
end

function y1 = MS(x,y,h,i)
    temp0 = f(x(i-2),y(i-2));
    temp1 = f(x(i-1),y(i-1));
    y1 = (3/(3+4*h)) * (y(i-2)+(h/3) * ( x(i)*exp( -4*x(i) ) +4*temp1+temp0));
end

function y1 = RK2(y0,h,x)
    k1 = f(x,y0);
    k2 = f(x+h/2,y0+h*k1/2);
    y1 = y0 + h*k2;
end

function y1 = f(x0,y0)
    y1 = x0*exp(-4*x0)-4*y0;
end