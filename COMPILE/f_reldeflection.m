function [ X_posinflpoint,reldef_mag,reldef_mag_sag,reldef_mag_hog,L_sag,L_hog] =...
    f_reldeflection(controller_int,X,Y,Uz)
%F_RELDEFLECTION This function compute the relative deflection of a curve 
% X,Y vector containing x,y coordinates of the pile heads
% Ux,z vector containing displacement of the pile heads in x,y direction

%% Interpolation with smoothingspline of Ux and Uz
% % It is possible to deactivate this part if controller_int is ==0

if controller_int==1
[xDataz, yDataz] = prepareCurveData( X, Uz);
% Set up fittype and options.
ft = fittype( 'smoothingspline' );
opts = fitoptions( ft );
opts.SmoothingParam = 0.5;

% Fit model to data.
[fitresultz, gofz] = fit( xDataz, yDataz, ft, opts );
% Plot fit with data.
figure(102);hold on
plot( fitresultz, xDataz, yDataz );

X_temp=min(X):(X(2)-X(1))/25:max(X);
Y_temp=X_temp*0;
X=X_temp;
Y=Y_temp;

Uz=fitresultz(X)';
end
%%

% Inflection point positions = position change of sign of Uz second derivative  
iider_fr=gradient(gradient(Uz)); %Uz second derivative 
cross = find((abs(diff(sign(iider_fr)))== 2) | (abs(diff(sign(iider_fr)))== 1)) ;% sign change and zero crossings
for ii=1:length(cross)% posinflpoint index of inflection points respect X Y vectors
    if abs(iider_fr(cross(ii)))<abs(iider_fr(cross(ii)+1))
        posinflpoint(ii)=cross(ii);
    else
        posinflpoint(ii)=cross(ii)+1;
    end
end
if exist('posinflpoint')==1 
    if and(min(posinflpoint)==1,max(posinflpoint)==length(Uz))
        posinflpoint=[posinflpoint];
    elseif max(posinflpoint)==length(Uz)
        posinflpoint=[1,posinflpoint];
    elseif min(posinflpoint)==1
        posinflpoint=[posinflpoint,length(Uz)];
    else
        posinflpoint=[1,posinflpoint,length(Uz)];
    end
else
posinflpoint=[1,length(Uz)];
end

X_posinflpoint=X(posinflpoint); %x coordinates inflection points

U_def_X=X';%coordinates deformed curve
U_def_Y=Y'+Uz';%coordinates deformed curve
mapxy=[U_def_X U_def_Y];%coordinates deformed curve
curvexy=mapxy(posinflpoint,:);%coordinates inflection points

[xy,distance,t_a] =distance2curve(curvexy,mapxy,'linear');%function to calculate relative deflection of each pile each (distance)

for i=1:length(posinflpoint)-1
    zone_sign(i)=sum(iider_fr(posinflpoint(i):posinflpoint(i+1))); %zone_sign component if + (-) indicates sagging (hogging)
    if zone_sign(i)<0   %deflection negative in sagging zone
        deflection(posinflpoint(i):posinflpoint(i+1))=-distance(posinflpoint(i):posinflpoint(i+1));
        deflection_mag(i)=min(deflection(posinflpoint(i):posinflpoint(i+1)));
    else                %deflection positive in hogging zone
        deflection(posinflpoint(i):posinflpoint(i+1))=distance(posinflpoint(i):posinflpoint(i+1));
        deflection_mag(i)=max(deflection(posinflpoint(i):posinflpoint(i+1)));
    end
end

for i=1:size(curvexy,1)-1 %lengthZone vector containing length of hogging and sagging zones
   lengthZone(i)=pdist(curvexy(i:i+1,:),'euclidean'); 
end



reldef_mag=deflection_mag./lengthZone;
if max(reldef_mag)>0
    reldef_mag_hog=max(reldef_mag);
else
    reldef_mag_hog=0;
end
if min(reldef_mag)<0
    reldef_mag_sag=min(reldef_mag);
else
    reldef_mag_sag=0;
end

[r,i_hog] = find(reldef_mag==max(reldef_mag));
[r,i_sag] = find(reldef_mag==min(reldef_mag));
if numel(i_hog)>0
    L_hog=lengthZone(i_hog(1));
else
    L_hog=NaN(1);
end
if numel(i_sag)>0
    L_sag=lengthZone(i_sag(1));
else
    L_sag=NaN(1);
end

figure(103); hold on
h(1)=plot(mapxy(:,1),mapxy(:,2),'r*');
h(2)=plot(curvexy(:,1),curvexy(:,2),'k-o');
h(3)=plot(X,deflection,'m--');
h(4)=plot(xy(:,1),xy(:,2),'g*');
h(5)=plot([min(X),max(X)],100*[reldef_mag_hog,reldef_mag_hog],'b-');
h(6)=plot([min(X),max(X)],100*[reldef_mag_sag,reldef_mag_sag],'b-');
line([mapxy(:,1),xy(:,1)]',[mapxy(:,2),xy(:,2)]','color',[0 0 1]);
set(gca,'FontName','Times New Roman','FontSize',9,'linewidth',1.0,'YDir','reverse','TickDir','out');
legend (h([3,6]),'Deflection','Relative Deflection*100','Location','Southeast')
hold off
end

