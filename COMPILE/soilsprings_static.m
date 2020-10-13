function [ kh ,kv ] = soilsprings_static(rm, Es,Ep,d,ni )

% Rigidezza delle molle del terreno secondo Dobry Gazetas alla profondità
% z=zgauss(ii)


% kh=1.2*Es;
% ch=2*d*ro*Vs*(1+3.4/(pi*(1-ni)))*omega^(.75)*(d/Vs)^(-.25)+eta*1.2*Es; 
% kv=0.6*Es*(1+.5*(omega*d/Vs)^.5);
% cv=1.2*omega^(.75)*(d/Vs)^(-.25)*pi*d*ro*Vs+eta*0.6*Es*(1+.5*(omega*d/Vs)^.5);
Gs=Es/2/(1+ni);
Ip=pi*(d/2)^4/4;

kh=0.65*Es/d/(1-ni^2)*(d^4*Es/Ep/Ip)^(1/12);
ch=0; 
kv=2*pi*Gs/log(rm/(d/2));
cv=0;


end

