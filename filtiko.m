function [ s ] = filtiko(L, ngrid,gradruido)

imag=i;
alfa=sqrt(gradruido);
delt=L/ngrid;
y=linspace(log(delt), log(L), ngrid);
%y=(1/L).*linspace(1,ngrid,ngrid);
s=zeros(1,ngrid);

for is=1:ngrid
    x=y(is);
    %terexp=2*power(2,2pi*imag*x);
    pote=2*pi*imag*x;
    terexp=2*power(2,pote);
    s(is)=(1/(2*terexp-1))*(1/(1+alfa*abs(x)));
    
end



end