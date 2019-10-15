function [a_xx,a_xy,a_yx,a_yy]=genDynamicCoe1(m,N,phi_st,phi_ex,Kr)

% tt=0:dt:T;
mm=1:m+1;
phi=zeros(N,length(mm));
for ii=1:N
    phi(ii,:)=mm*2*pi/m/N+(ii-1)*2*pi/N;
end
g=zeros(N,length(mm));
for ii=1:N
    for jj=1:length(mm)
        if phi(ii,jj)>=phi_st && phi(ii,jj)<=phi_ex
            g(ii,jj)=1;
        end
    end
end

a_xx=sum(-g.*(sin(2*phi)+Kr*(1-cos(2*phi))));
a_xy=sum(-g.*((1+cos(2*phi))+Kr*sin(2*phi)));
a_yx=sum(g.*((1-cos(2*phi))-Kr*sin(2*phi)));
a_yy=sum(g.*(sin(2*phi)-Kr*(1+cos(2*phi))));
