%%
%最小二乘法 （一）
x = [2.1 2.5 2.8 3.2];
y = [0.6087 0.6849 0.7368 0.8111];
m = 4;
for i=1:m
   y_1(i) = y(i)*x(i);
   y_2(i) = y(i)*x(i)^2;
   x_1(i) = y(i)^2;
   x_2(i) = y(i)^2*x(i);
   x_3(i) = y(i)^2*x(i)^2;
end
A = [sum(x_1) sum(x_2);sum(x_2) sum(x_3)];
Y = [sum(y_1);sum(y_2)];
X = A\Y;
disp(X)
for i=1:m
   y_sol(i) = x(i)/(X(1)+X(2)*x(i));
end
if  X(2)>= 0
    fprintf("function is x/(%.16f+%.16fx)\n",X(1),X(2));
else
    fprintf("function is x/(%.16f%.16fx)\n",X(1),X(2));
end
fprintf("Error-2-Norm is %.16f\n",norm(y_sol-y,2));
%%
%最小二乘法 （二）
x = [2.1 2.5 2.8 3.2];
y = [0.6087 0.6849 0.7368 0.8111];
m = 4;
for i=1:m
   y_1(i) = y(i)*x(i);
   y_2(i) = y(i)^2*x(i);
   x_1(i) = x(i)^2;
   x_2(i) = y(i)*x(i)^2;
   x_3(i) = y(i)^2*x(i)^2;
end
A = [sum(x_1) -sum(x_2);sum(x_2) -sum(x_3)];
Y = [sum(y_1);sum(y_2)];
alpha = A\Y;
X=[1/alpha(1);alpha(2)/alpha(1)];
for i=1:m
   y_sol(i) = x(i)/(X(1)+X(2)*x(i));
end
if  X(2)>= 0
    fprintf("function is x/(%.16f+%.16fx)\n",X(1),X(2));
else
    fprintf("function is x/(%.16f%.16fx)\n",X(1),X(2));
end
fprintf("Error-2-Norm is %.16f\n",norm(y_sol-y,2));

