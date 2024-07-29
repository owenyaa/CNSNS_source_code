function [Xmin] = getmin(y)
a=length(y);
j=1;
for i=(a-1)/2:a
    b=(y(i,1)-y(i-2,1))/2;
    c=(y(i,1)+y(i-2,1))/2-y(i-1,1);
    if y(i-2,1)>=y(i-1,1)&&y(i-1,1)<=y(i,1)&&c==0
        Xmin(j)=y(i-1,1);
        j=j+1;
    elseif y(i-2,1)>=y(i-1,1)&&y(i-1,1)<=y(i,1)
        Xmin(j)=y(i-1,1)+b^2/(4*c);
        j=j+1;
    end
end