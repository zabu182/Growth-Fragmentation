function [ y ] = B( v )

stemp=size(v);
s=stemp(2);
y=v;
for i=1:s
    x=v(i);
    y(i)=(x^2)*exp(-(x+1/x));
    %y(i)=exp(-(x-10)*(x-10));
    

    
end




end