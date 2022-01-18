function [ y ] = truncat( v,ordertru )
stemp=size(v);
s=stemp(2)/(power(2,ordertru));
y=zeros(1,s);


for i=1:s
    
    y(i)=v(i);
    
end

end