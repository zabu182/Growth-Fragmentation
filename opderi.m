function [ y ] = opderi( v,dife )
stemp=size(v);
s=stemp(2);
y=zeros(stemp);

%y(s)=(v(s)-v(s-1))/dife;
y(1)=(v(3)-4*v(2)+3*v(1))/(-2*dife);
y(2)=(v(4)-4*v(3)+3*v(2))/(-2*dife);
y(s)=(2*v(s-2)-v(s-1)-v(1))/-3*dife;
y(s-1)=(2*v(s-3)-v(s-2)-2*v(s-1))/-3*dife;


for i=3:s-2
    
    y(i)=(-v(i+2)+8*v(i+1)-8*v(i-1)+v(i-2))/(12*dife);
    
end

end