function [ ux,uz ]=u_VER(z,x,h,R,epsilon,delta,ni)  

%x,z are the spatial coordinates
%R is the tunnel radius
%h is the tunnel depth
%epsilon is equal to Vlt./2
%ni Poisson's ratio


z1=z-h;
z2=z+h;
m=1./(1-2.*ni);
k=ni./(1-ni); 
r1=(x.^2+z1.^2).^.5;
r2=(x.^2+z2.^2).^.5;

u_x_ue1=-epsilon.*R.^2.*(x./r1.^2+x./r2.^2);
u_x_ud1=+delta.*R.^2.*(x.*(x.^2-k.*z1.^2)./r1.^4+x.*(x.^2-k.*z2.^2)./r2.^4);
u_y_ue1=-epsilon.*R.^2.*(z1./r1.^2+z2./r2.^2);
u_y_ud1=+delta.*R.^2.*(z1.*(k.*x.^2-z1.^2)./r1.^4+z2.*(k.*x.^2-z2.^2)./r2.^4);

u_x_ue2=-2.*epsilon.*R.^2.*x./m.*(1./r2.^2-2.*m.*z.*z2./r2.^4);
u_x_ud2=-4.*delta.*R.^2.*x.*h./(m+1).*(z2./r2.^4+m.*z.*(x.^2-3.*z2.^2)./r2.^6);
u_y_ue2=2.*epsilon.*R.^2./m.*((m+1).*z2./r2.^2-m.*z.*(x.^2-z2.^2)./r2.^4);
u_y_ud2=-2.*delta.*R.^2.*h.*((x.^2-z2.^2)./r2.^4+m./(m+1).*2.*z.*z2.*(3.*x.^2-z2.^2)./r2.^6);

ux=u_x_ue1+u_x_ue2+u_x_ud1+u_x_ud2;
uz=u_y_ue1+u_y_ue2+u_y_ud1+u_y_ud2;


