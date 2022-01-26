clear;
real = 7.9549265210128452745;
for i=1:15
    err(i) = Composite(i,real);
end
figure
semilogy([1:15],err);

function err = Composite(M,real)
    h = 2*pi/M;
    sum = 0;
    x(1) = -pi;
    for i=2:M
        x(i) = x(i-1)+h;
        sum = sum + f(x(i)) ;
    end
    result = h*(0.5*(f(-pi))+sum+0.5*(f(pi)));
    err = abs(result - real);
end
function y=f(x)
    y=exp(cos(x));
end