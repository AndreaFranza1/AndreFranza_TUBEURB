function [ eps_ht,eps_hc,pos_ht,pos_hc] =...
    f_hrzstr(controller_int,X,Y,Ux)
%F_RELDEFLECTION This function compute the relative deflection of a curve 
% X,Y vector containing x,y coordinates of the pile heads
% Ux,z vector containing displacement of the pile heads in x,y direction

%% Interpolation with smoothingspline of Ux and Uz
% % It is possible to deactivate this part

if controller_int==1
[xDatax, yDatax] = prepareCurveData( X, Ux);
% Set up fittype and options.
ft = fittype( 'smoothingspline' );
opts = fitoptions( ft );
opts.SmoothingParam = 0.5;

% Fit model to data.
[fitresultx, gofx] = fit( xDatax, yDatax, ft, opts );
% Plot fit with data.
figure(101);hold on
plot( fitresultx, xDatax, yDatax );

X_temp=min(X):(X(2)-X(1))/25:max(X);
Y_temp=X_temp*0;
X=X_temp;
Y=Y_temp;

Ux=fitresultx(X)';
end
%%
eps=gradient(Ux,X(2)-X(1));
[eps_ht,pos_ht]=(max(eps));
[eps_hc,pos_hc]=(min(eps));

X_pos_eps_ht=X(pos_ht);
X_pos_eps_hc=X(pos_hc);


if max(eps)>0
    eps_ht=max(eps);
else
    eps_ht=0;
end
if min(eps)<0
    eps_hc=min(eps);
else
    eps_hc=0;
end


figure(104); hold on

h(2)=plot(X,eps,'m--');
h(1)=plot(X,Ux/10,'r*')
h(3)=plot([X_pos_eps_ht,X_pos_eps_ht],[-eps_ht,eps_ht],'k-');
h(4)=plot([X_pos_eps_hc,X_pos_eps_hc],[-eps_hc,eps_hc],'k-');
h(5)=plot([min(X),max(X)],[eps_ht,eps_ht],'b-');
h(6)=plot([min(X),max(X)],[eps_hc,eps_hc],'b-');
set(gca,'FontName','Times New Roman','FontSize',9,'linewidth',1.0,'YDir','reverse','TickDir','out');
legend (h([1,2]),'Horizontal Movements/10','Horizontal strains','Southeast')
hold off


end

