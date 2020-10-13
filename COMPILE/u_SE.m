function [ ux,uz ]=u_SE(z,x,h,R,epsilon,id,reg_ca,reg_cb,reg_c1,reg_c2,reg_c3,reg_c4,reg_c5,reg_c6,reg_cax,reg_cbx,reg_c1x,reg_c2x,reg_c3x,reg_c4x,reg_c5x,reg_c6x)  

%x,z are the spatial coordinates
%R is the tunnel radius
%h is the tunnel depth
%epsilon is equal to Vlt./2
%alfa is equal to 1 for elastic behaviour
%delta is equalt to epsilon hypotesis

%Calibration assumption


%Parameter
zt=h;
cd=(zt-R)./2./R;
vlt_p=epsilon.*2.*100;
beta=h./R;% beta=3;
ro=1;

%Model dimensionaless coordinates
xp=x./h;
z1 = z-h;
z2 = z+h;
zp1=z1./h;
zp2=z2./h;
zp=z./h;

rp1=((x.^2+z1.^2).^.5)./h;
rp2=((x.^2+z2.^2).^.5)./h;

r1=((x.^2+z1.^2).^.5);
r2=((x.^2+z2.^2).^.5);


%% Zita definition

% CA = ma.*vlt_p+ qa;
% CB = mb.*vlt_p+ qb;
% c1 = m1.*vlt_p+ q1;
% c1x= m1x.*vlt_p + q1x;
% c2 = m2.*vlt_p+ q2 ;
% c3 = m3.*vlt_p+ q3;
% c4 = m4.*vlt_p+ q4;
% c5 = m5.*vlt_p+ q5;

CA = reg_ca(1).*vlt_p+ reg_ca(2);
CB = reg_cb(1).*vlt_p+ reg_cb(2);
c1 = reg_c1(1).*vlt_p+ reg_c1(2);
c2 = reg_c2(1).*vlt_p+ reg_c2(2) ;
c3 = reg_c3(1).*vlt_p+ reg_c3(2);
c4 = reg_c4(1).*vlt_p+ reg_c4(2); 
c5 = reg_c5(1).*vlt_p+ reg_c5(2);
c6 = reg_c6(1).*vlt_p+ reg_c6(2) ;

CAx = reg_cax(1).*vlt_p+ reg_cax(2);
CBx = reg_cbx(1).*vlt_p+ reg_cbx(2);
c1x = reg_c1x(1).*vlt_p+ reg_c1x(2);
c2x = reg_c2x(1).*vlt_p+ reg_c2x(2) ;
c3x = reg_c3x(1).*vlt_p+ reg_c3x(2);
c4x = reg_c4x(1).*vlt_p+ reg_c4x(2); 
c5x = reg_c5x(1).*vlt_p+ reg_c5x(2);
c6x = reg_c6x(1).*vlt_p+ reg_c6x(2) ;

zitaz = CA.*exp(-(c1.*zp.^2+c2.*xp.^2+c6.*xp.^4))...
        +CB.*exp(-(c3.*(zp-c4).^2+c5.*xp.^2));
zitax = CAx.*exp(-(c1x.*zp.^2+c2x.*xp.^2+c6x.*xp.^4))...
        +CBx.*exp(-(c3x.*(zp-c4x).^2+c5x.*xp.^2));
    
    
% % 
% ma=(0.04763.*id-0.15)..*exp(-(cd-1)..^6.05./100);
% qa=(-0.78.*id+1.18)..*exp(-(cd-1)..^6.5./1000)+0.3;
% mb=(0.18)..*exp(-(cd-1)..^7.48./15);
% qb=0;
% m1=(-0.5..*id+0.66).*exp(-(cd-1)..^2.3./4)+0.09; %
% q1=cd../cd.*(0..*id+0)+1.0;
% m2=(-0.1..*id+0.263).*exp(-(cd-1)..^3.957./25);
% q2=0;
% m3=(1.244..*id+0.016).*exp(-(cd-1)..^3.95./6.16)+0.0; %
% q3=0;
% m4=0;
% q4=(cd)../(cd+0.5);
% m5=exp(-(cd-1)..^3.582./25).*exp(2.507.*id)+0;
% q5=0;
% m1x=0;
% q1x=log((cd)..^5.7)+3.3+id..*0;
% % 






%% Elastic and semianalitical ground movements

ux_elast= 2.*epsilon.*R.*(1./beta).*   (-((xp.*(1 - (ro.*(xp.^2 - zp1.^2))./rp1.^2))./(2.*rp1.^(2))) -     (xp.*(1 - (ro.*(xp.^2 - zp2.^2))./rp2.^2))./(2.*rp2.^(2)) + ...
    (4.*xp.*zp.*(zp2./rp2.^2 - (ro.*(xp.^2 - 3.*zp2.^2))./rp2.^4))./     (2.*rp2.^(2)));

uz_elast= 2.*epsilon.*R.*(1./beta).* (-((zp1.*(1 - (ro.*(xp.^2 - zp1.^2))./rp1.^2))./(2.*rp1.^(2))) +(zp2.*(1 + (ro.*(xp.^2 - zp2.^2))./rp2.^2))./(2.*rp2.^(2)) - ...
    ((2.*(zp + ro).*(xp.^2 - zp2.^2))./rp2.^2 +     (4.*ro.*zp.*zp2.*(3.*xp.^2 - zp2.^2))./rp2.^4)./(2.*rp2.^(2)));



%General expression for incompressible elastic soil Gonzales
% ux_elast= 2.*epsilon.*R.*(R./h).^(2.*alfa - 1).*   (-((xp.*(1 - (ro.*(xp.^2 - zp1.^2))./rp1.^2))./(2.*rp1.^(2.*alfa))) -     (xp.*(1 - (ro.*(xp.^2 - zp2.^2))./rp2.^2))./(2.*rp2.^(2.*alfa)) + ...
%     (4.*xp.*zp.*(zp2./rp2.^2 - (ro.*(xp.^2 - 3.*zp2.^2))./rp2.^4))./     (2.*rp2.^(2.*alfa)));
% uz_elast = 2.*epsilon.*R.*(R./h).^(2.*alfa - 1).* (-((zp1.*(1 - (ro.*(xp.^2 - zp1.^2))./rp1.^2))./(2.*rp1.^(2.*alfa))) +(zp2.*(1 + (ro.*(xp.^2 - zp2.^2))./rp2.^2))./(2.*rp2.^(2.*alfa)) - ...
%     ((2.*(zp + ro).*(xp.^2 - zp2.^2))./rp2.^2 +     (4.*ro.*zp.*zp2.*(3.*xp.^2 - zp2.^2))./rp2.^4)./(2.*rp2.^(2.*alfa)));

% ux_elast=-1.*2.*epsilon.*R.^2.*(x./2./r1.^2.*(1-(x.^2-z1.^2)./r1.^2)+x./2./r2.^2.*(1-(x.^2-z2.^2)./r2.^2)-4.*x.*z./2./r2.^4.*(z2-h.*(x.^2-3.*z2.^2)./r2.^2));
%  
% uz_elast=-1.*2.*epsilon.*R.^2.*(z1./2./r1.^2.*(1-(x.^2-z1.^2)./r1.^2)-z2./2./r2.^2.*(1+(x.^2-z2.^2)./r2.^2)+1./2./r2.^4.*(2.*(z+h).*(x.^2-z2.^2)+4.*h.*z.*z2.*(3.*x.^2-z2.^2)./r2.^2));

% ux_elast = 2.*epsilon.*R.^2.*((-(x./(2.*r1.^2))).*(1 - (x.^2 - z1.^2)./r1.^2) -(x./(2.*r2.^2)).*(1 - (x.^2 - z2.^2)./r2.^2) + ((4.*x.*z)./(2.*r2.^4)).* (z2 - (h.*(x.^2 - 3.*z2.^2))./r2.^2));
% 
% uz_elast = 2.*epsilon.*R.^2.*((-(z1./(2.*r1.^2))).*(1 - (x.^2 - z1.^2)./r1.^2) + (z2./(2.*r2.^2)).*(1 + (x.^2 - z2.^2)./ r2.^2) - ((2.*(z + h).*(x.^2 - z2.^2))./(2.*r2.^4) + (4.*h.*z.*z2.*(3.*x.^2 - z2.^2))./(2.*r2.^6)));


ux = zitax .*  ux_elast;

uz = zitaz .*  uz_elast;


