%%% Inverse problem Growth Fragmentation

%% Main parameters

clc
clear all
close all



ngrid=500; %% Ngrid variable x
L=20; %lenght interval variable x
delt1=L/ngrid;
ngridtime=10000; %% Ngrid variable t
L2=200; %lenght interval variable t
delt2=L2/ngridtime;
ontes=ones(1,ngrid);
basx=delt1.*linspace(1,ngrid,ngrid);


%----------End parameter-------------

MatB=B(basx);
%MatB=2*ones(1,power(2,ngridtime)*ngrid);
Matg=g(basx);
%Matg=ones(1,power(2,ngridtime)*ngrid);
n=n0(basx);
%n=ones(1,power(2,ngridtime)*ngrid);

%--------Data set--------------------

%%-------------- Direct Problem-------

for i=2:ngridtime+1
    
    
    Fpart=truncat( Matg,0 ).*truncat( n,0);
    Fpart=opderi(Fpart,delt1);
    Spart=truncat( MatB,0 ) .*truncat( n,0);
    Tpart=opkernel( MatB,n );
    
    n=truncat( n,0)+ delt2.*(Tpart-Fpart-Spart);
    n=n*(1/(delt1*norm(n,1)));

    
end
nsol=n;
nsolcop=n;
Hsol=MatB.*n;
Hsolco=Hsol;

eigval=dot(ontes,(Tpart-Fpart-Spart))/(dot(ontes,n));
eigvalcop=eigval;
    


%%%%---------End Direct Problem

epsima=linspace(power(10,-4),0.5,200);
eti=zeros(1,200);
elan=eti;
equa=eti;
eti2=eti;
elan2=eti;
equa2=eti;
for coepsi=1:200
   
    gradruido=epsima(1,coepsi);
    close all
    clc
    n=nsolcop;
    nsol=n;
    Hsol=Hsolco;
    eigval=eigvalcop;


%%%%%%Inverse Problem%%%%%%%%%%%%%%%%%%
yno=rand(1,ngrid);
yno=yno/(delt1*norm(yno,2));



%%%% Noise %%%%
%n=n+gradruido*yno.*n;
n=n+gradruido*yno
n=(abs(n)+n)/2;
eigval=eigval+gradruido;%.*yno1;

%%%% end noisse%%%%%%
Fpart=truncat( Matg,0 ).*truncat( n,0);
Fpart=opderi(Fpart,delt1);


Invpriter1=filtiko(L, ngrid,gradruido).*newfft(Fpart,L, ngrid);

Invpriter2=filtland(L, ngrid,gradruido).*newfft(Fpart,L, ngrid);

Invsegter=filtstand(L, ngrid,gradruido).*newfft(n,L, ngrid);

Invpriter31=filtquasere(L, ngrid,gradruido).*newfft(Fpart,L, ngrid);

Invsegter32=filtquasere(L, ngrid,gradruido).*newfft(n,L, ngrid);

H1sol=newifft(Invpriter1+eigval.*Invsegter,L, ngrid);
H2sol=newifft(Invpriter2+eigval.*Invsegter,L, ngrid);
H3sol=newifft(Invpriter31+eigval.*Invsegter32,L, ngrid);

%%%---------------End Inverse problem-----------

%%%--------------Plot figures---n,H
x=truncat(basx,1);
nsol=truncat(nsol,1);
H1sol=real(truncat(H1sol,1));
H2sol=real(truncat(H2sol,1));
H3sol=real(truncat(H3sol,1));
Hsol=truncat(Hsol,1);

xa=x;
x=x*(3/10);
na=nsol;
H1a=H1sol;
H2a=H2sol;
H3a=H3sol;
Ha=Hsol;

nsol=interp1(xa,nsol,x);
H1a=interp1(xa,H1a,x);
H2a=interp1(xa,H2a,x);
H3a=interp1(xa,H3a,x);
Ha=interp1(xa,Ha,x);


figure
plot(x,nsol,'k')
title('N(x)')
hold on

 
figure
plot(x,Ha,'k',x,H1a,'b--',x,H2a,'m--',x,H3a,'r--')
h=legend('Exact numerical solution','Tikhonov filtering','Landweber filtering','Quasi Reversibility Method','location', 'northeast' )
title('H(x)')
set(h(1),'linewidth',2);
hold on 

x=xa;
nsol=na;



%%-------fin


%%---------Plot B



H1sol(isnan(H1sol))=0;
H2sol(isnan(H2sol))=0;
H3sol(isnan(H3sol))=0;


Btol=H1sol./nsol;
Bto2=H2sol./nsol;
Bto3=H3sol./nsol;



Btol=(Btol+abs(Btol))/2;
Bto2=(Bto2+abs(Bto2))/2;
Bto3=(Bto3+abs(Bto3))/2;

Btol= medfilt1(Btol,30);
Bto2= medfilt1(Bto2,30);
Bto3= medfilt1(Bto3,30);
Bsol=truncat(MatB,1);
xa=x;
x=x*(3/10);

Bsol=interp1(xa,Bsol,x);
Btol=interp1(xa,Btol,x);
Bto2=interp1(xa,Bto2,x);
Bto3=interp1(xa,Bto3,x);

%Btol=truncat(Btol,2);
%Bto2=truncat(Bto2,2);
%Bto3=truncat(Bto3,2);

%x=truncat(x,2);

Bsol(isnan(Bsol))=0;
Btol(isnan(Btol))=0;
Bto2(isnan(Bto2))=0;
Bto3(isnan(Bto3))=0;



eti(1,coepsi)=norm(Hsol-H1sol)
elan(1,coepsi)=norm(Hsol-H2sol)
equa(1,coepsi)=norm(Hsol-H3sol)

eti2(1,coepsi)=norm(Bsol-Btol)
elan2(1,coepsi)=norm(Bsol-Bto2)
equa2(1,coepsi)=norm(Bsol-Bto3)

norm(Bsol-Btol)
norm(Bsol-Bto2)
norm(Bsol-Bto3)


figure
h=plot(x,Bsol,'k',x,Btol,'b--',x,Bto2,'m--',x,Bto3,'r--')
legend('Exact solution','Tikhonov filtering','Landweber filtering','Quasi Reversibility Method','location', 'Northwest' )
title('B(x)')
set(h(1),'linewidth',2);
hold on 

end
figure
h=plot(log(epsima),log(delt1*eti),log(epsima),log(delt1*elan), log(epsima),log(delt1*equa))
legend('Tikhonov filtering','Landweber filtering','Quasi Reversibility Method','location', 'Northwest' )
title('Numerical error for the function H')
xlabel('Log(\epsilon)') 
ylabel('L^2 error in logarithmic scale') 
set(h(1),'linewidth',2);
set(h(2),'linewidth',2);
set(h(3),'linewidth',2);

figure
h=plot(log(epsima),log(delt1*eti2),log(epsima),log(delt1*elan2), log(epsima),log(delt1*equa2))
legend('Tikhonov filtering','Landweber filtering','Quasi Reversibility Method','location', 'Northwest' )
title('Numerical error for the function B')
xlabel('Log(\epsilon)') 
ylabel('L^2 error in logarithmic scale') 
set(h(1),'linewidth',2);
set(h(2),'linewidth',2);
set(h(3),'linewidth',2);