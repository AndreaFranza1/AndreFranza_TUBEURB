function [ uxij_3,uyij_3,uzij_3 ] = int_factor_mindlin_j_3_vect( Xi,Yi,Zi,Xj,Yj,Zj,Gs,nus,kh,kv,d )
% uxij_3,uyij_3,uzij_3 displacement at point i due to unit force in j
% in direction 3 (vertical one)
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



t1w=(3-4.*nus)./R1;
t2w=(8.*(1-nus).^2-(3-4.*nus))./R2;
t3w=(c-z).^2./R1.^3;
t4w=((3-4.*nus).*(z+c).^2-2.*z.*c)./R2.^3;
t5w=6.*c.*z.*(z+c).^2./R2.^5;
w=1./16./pi./Gs./(1-nus).*(t1w+t2w+t3w+t4w+t5w);


t1U=(z-c)./R1.^3;
t2U=(3-4.*nus).*(z-c)./R2.^3;
t3U=+6.*c.*z.*(z+c)./R2.^5;  
t4U=-4.*(1-nus).*(1-2.*nus)./R2./(R2+z+c); 
U=1./16./pi./Gs./(1-nus).*(t1U+t2U+t3U+t4U);

w(isnan(w(:,1)))=1/kv;



uxij_3=U.*x;%Mindlin paper
uyij_3=U.*y;%Mindlin paper

uxij_3(isnan(uxij_3(:,1)))=0;
uyij_3(isnan(uyij_3(:,1)))=0;

uzij_3=w;


end

