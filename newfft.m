function [ tra] = newfft(n,L, ngrid)

imag=i;
y=zeros(1,ngrid);
delt=L/ngrid;
malha=linspace(delt, L, ngrid);
malhafrecueuniespa=linspace(log(delt),log(L),ngrid);
delfre=(log(L)-log(delt))/ngrid;




for ik=1:ngrid
    


    y(ik)=log(malha(ik));
    
end



vq = interp1(y,n,malhafrecueuniespa);
tra=zeros(1,ngrid);


for ik=1:ngrid
    
xi=malhafrecueuniespa(1,ik);
sum=0;
for ik2=1:ngrid
    x=malhafrecueuniespa(1,ik2);


    sum=sum+vq(ik2)*exp(-2*pi*imag*x*xi);
    
end

    tra(1,ik)=sum*delfre;
end





end