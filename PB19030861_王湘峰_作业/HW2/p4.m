syms a b;
x0=[2.1 2.5 2.8 3.2];
y0=[0.6087 0.6849 0.7368 0.8111];
eq1=0;eq2=0;
for i=1:4
    eq1=eq1+(x0(i)-y0(i)*(a+b*x0(i)))*y0(i);
    eq2=eq2+(x0(i)-y0(i)*(a+b*x0(i)))*y0(i)*x0(i);
end
s=solve(eq1,eq2,a,b);
fprintf('a=%.14f\nb=%.14f\n',s.a,s.b);
%计算误差2-范数
err=0;
for i=1:4
    err=err+(y0(i)-x0(i)/(s.a+s.b*x0(i)))^2;
end
err=err^0.5;
fprintf('误差的2-范数为%.14f\n',err);