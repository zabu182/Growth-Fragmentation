function [ s ] = filtland(L, ngrid,gradruido)

imag=i;
m=10;
base=(2*power(m,m+1)*(1/gradruido));
alfa=power(base,(2/((2*m)+1)));
delt=L/ngrid;
y=linspace(log(delt), log(L), ngrid);
%y=(1/L).*linspace(1,ngrid,ngrid);
s=zeros(1,ngrid);

for is=1:ngrid
    x=y(is);
    %terexp=2*power(2,2pi*imag*x);
    pote=2*pi*imag*x;
    terexp=2*power(2,pote);
    inteter=1-(1/(1+power(x,2)));
    s(is)=(1/(2*terexp-1))*(1-power(inteter,alfa));
    
end



end