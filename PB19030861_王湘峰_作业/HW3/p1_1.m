H=zeros(1,16);
y=zeros(1,16);
exact=cos(1.2);
err=zeros(1,16);
for i=0:15
    H(i+1)=10^(-i);
    y(i+1)=forward(1.2,H(i+1));
end
for i=1:16
    err(i)=abs(y(i)-exact);
end
loglog(H,err);


function y=f(x)
    y=sin(x);
end
function y=forward(x,h)
    y=(f(x+h)-f(x))/h;
end