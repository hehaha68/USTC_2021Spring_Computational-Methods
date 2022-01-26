clear;
real = cos(1.2);
n = [1:40];
N = Table(1.2,67,40);
for i=1:40
    err(i) = abs(N(i,i)-real);
end
figure
semilogy(n,err);
fprintf("导数值为：%.16f。\n",N(14,14));
fprintf("误差为：%.16g。\n",abs(N(14,14)-real));
fprintf("h初始值为：67。\n");
fprintf("外推次数为：14。\n");

function N = Table(x0,h,n)
    N = zeros(n,n);
    for i = 1:n
        N(i,1) = (sin(x0+h/2^(i-1))-sin(x0))/(h/2^(i-1));
    end
    for i = 2:n
        for j = i:n
            N(j,i) = N(j,i-1) + (N(j,i-1) - N(j-1,i-1))/(2^(i-1)-1);
        end
    end
end

