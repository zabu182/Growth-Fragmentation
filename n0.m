function [ y ] = n0( v )

stemp=size(v);
s=stemp(2);

for i=1:s
    x=v(i);
    y(i)=power(x-8,2)*exp(-(x-8)*(x-8)/2);

    

    
end



end