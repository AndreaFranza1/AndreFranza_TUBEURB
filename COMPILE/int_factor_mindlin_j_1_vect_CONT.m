function [ uxij_1,uyij_1,uzij_1 ] = int_factor_mindlin_j_1_vect_CONT( Xi,Yi,Zi,Xj,Yj,Zj,Gs,nus )
% uxij_1,uyij_1,uzij_1 displacement at point i due to unit force in j
% in direction 1 (x)
%Xi,Yi,Zi,Xj,Yj,Zj global coordinates of i and j point

z=Zi;
c=Zj;
x=Xi-Xj;
y=Yi-Yj;
r=((Xi-Xj).^2+(Yi-Yj).^2).^.5;
% 
% z=Zi_gauss;
% c=Zj_gauss;
% x=Xi_gauss-Xj_gauss;
% y=Yi_gauss-Yj_gauss;
% r=((Xi_gauss-Xj_gauss).^2+(Yi_gauss-Yj_gauss).^2).^.5;




R1=(r.^2+(c-z).^2).^.5;
R2=(r.^2+(c+z).^2).^.5;

t1u=(3-4.*nus)./R1;
t2u=1./R2;
t3u=x.^2./R1.^3;
t4u=(3-4.*nus).*x.^2./R2.^3;
t5u=2.*c.*z./R2.^3.*(1-3.*x.^2./R2.^2);
t6u=4.*(1-nus).*(1-2.*nus)./(R2+c+z).*(1-(x.^2./R2./(R2+c+z)));
u=1./16./pi./Gs./(1-nus).*(t1u+t2u+t3u+t4u+t5u+t6u);

t1v=1./R1.^3;
t2v=(3-4.*nus)./R2.^3;
t3v=-6.*c.*z./R2.^5;
t4v=-4.*(1-nus).*(1-2.*nus)./R2./(R2+z+c).^2;
v=x.*y./16./pi./Gs./(1-nus).*(t1v+t2v+t3v+t4v);

t1w=(z-c)./R1.^3;
t2w=(3-4.*nus).*(z-c)./R2.^3;
t3w=-6.*c.*z.*(z+c)./R2.^5;
t4w=+4.*(1-nus).*(1-2.*nus)./R2./(R2+z+c);
w=x./16./pi./Gs./(1-nus).*(t1w+t2w+t3w+t4w);


% u(isnan(u(:,1)))=1/kh/passo;
% v(isnan(v(:,1)))=0;
% w(isnan(w(:,1)))=0;

uxij_1=u;
uyij_1=v;
uzij_1=w;


end

