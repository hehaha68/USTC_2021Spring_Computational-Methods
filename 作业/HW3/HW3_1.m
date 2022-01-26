clear;
real = cos(1.2);
h = [1,10^(-1),10^(-2),10^(-3),10^(-4),10^(-5),10^(-6) ...
    ,10^(-7),10^(-8),10^(-9),10^(-10),10^(-11),10^(-12),10^(-13) ...
    ,10^(-14),10^(-15)];
for i=1:16
    diff = FD(1.2,h(i));
    fprintf("h为10^(%d)时，导数值为：%.16f\n",-(i-1),diff);
    err(i) = abs(real-diff);
end

figure
loglog(h,err);

function diff = FD(x0,h)
    diff = (sin(x0+h)-sin(x0))/h;
end