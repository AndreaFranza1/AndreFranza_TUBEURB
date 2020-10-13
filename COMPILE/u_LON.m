function [ ux,uz ]=u_LON(z,x,h,R,epsilon,delta,ni)  

%x,z are the spatial coordinates
%R is the tunnel radius
%h is the tunnel depth
%NB according to lonaghatan epsilon is epsilon_zero equal to
%epsilon*2=Vlt./2*2=Vlt
%ni Poisson's ratio

z1=z-h;
z2=z+h;
m=1./(1-2.*ni);
k=ni./(1-ni); 
r1=(x.^2+z1.^2).^.5;
r2=(x.^2+z2.^2).^.5;



uz=2.*exp(-1.*(1.38.*x.^2./(h+R).^2+0.69.*z.^2./h.^2)).*epsilon...
    .*R.^2.*(-1.*(z-h)./(x.^2+(z-h).^2)+(3-4.*ni).*(z+h)./(x.^2+(z+h).^2)...
    -(2.*z.*(x.^2-(z+h).^2)./(x.^2+(z+h).^2).^2));
ux=-1.*x.*2.*exp(-1.*(1.38.*x.^2./(h+R).^2+0.69.*z.^2./h.^2)).*epsilon...
    .*R.^2.*(1./(x.^2+(z-h).^2)+(3-4.*ni)./(x.^2+(z+h).^2)-(4.*z.*(z+h))./(x.^2+(z+h).^2).^2);


