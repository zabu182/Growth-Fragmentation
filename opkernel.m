function [ y ] = opkernel( C,n )
stemp=size(n);
s=stemp(2);
temp=floor(s/2);
y=zeros(1,s);


for i=1:temp
    
    y(i)=4*C(2*i)*n(2*i);
    
end

end