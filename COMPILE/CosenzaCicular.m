function [Ned,M_rd3] = CosenzaCicular(D,rho)

%NB N positive of compression

% % D=0.600;
c=0.045;
fcd= 14E6; %design value of concrete compressive strength
fyd=391E6; %design value of steel strength
% % rho=0.01;  %ratio steel/concrete


%%
x=1.389*D/5:D/1000:100*D;

R=D/2;
fpcd=0.9*fcd;
As=rho*pi*R^2;
omega=As/pi/R^2*fyd/fpcd;


teta=acos((R-0.8.*x)./x);

% figure(1)
% plot(x/R,teta/pi)


%%
Nc=R.^2/2*(2.*teta-sin(2.*teta)).*fpcd;
Ns=((teta/pi)-(1-teta/pi)).*As.*fyd;

Ned=Nc+Ns;
nu=Ned/pi/R.^2/fpcd;

Nmax=pi*R^2*fpcd+As.*fyd;
Nmin=-As.*fyd;

%in the paper equation there is an error 
% mu_rd=2/3/pi*sin(teta).^3+omega/pi.*(1-c/R).*sin(teta);
% M_rd=mu_rd*(2*pi*R^3*fpcd);
% M_rd2=4/3*R^3*sin(teta).^3*fpcd+2/pi*(R-c)*As.*sin(teta)*fyd;
% test1=M_rd./M_rd2;

F_c=R^2/2*(2.*teta-sin(2.*teta))*fpcd;
Br_c=4/3.*sin(teta).^3./(2.*teta-sin(2.*teta)).*R;
M_rd3=F_c.*Br_c+2/pi*(R-c)*As.*sin(teta)*fyd;

Ned=[Nmin,Ned,Nmax];
M_rd3=[0,M_rd3,0];

Ned=-Ned; %tensile axial forces are positive in the main script
end

