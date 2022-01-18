function [ y ] = g(v )

stemp=size(v);
s=stemp(2);
y=v;
for i=1:s
    x=v(i);
    y(i)=x*exp(-(x+1/x));


  
    

    
end



end