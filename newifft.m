function [ vq] = newifft(n,L, ngrid)

imag=i;
y=zeros(1,ngrid);
delt=L/ngrid;
malha=linspace(delt, L, ngrid);
malhafrecueuniespa=linspace(log(delt),log(L),ngrid);
delfre=(log(L)-log(delt))/ngrid;

%-------------------



re=zeros(1,ngrid);



for ik=1:ngrid
    
xi=malhafrecueuniespa(1,ik);
sum=0;
for ik2=1:ngrid
    x=malhafrecueuniespa(1,ik2);


    sum=sum+n(ik2)*exp(2*pi*imag*x*xi);
    
end

    re(1,ik)=(real(sum*delfre));
end
%-----------------


for iss=1:ngrid
    


    y(iss)=exp(malhafrecueuniespa(iss));
    
end




vq = interp1(y,re,malha);




end