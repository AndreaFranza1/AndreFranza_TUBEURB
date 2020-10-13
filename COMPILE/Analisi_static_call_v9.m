% This code run the analysis of tunnel pile structure interaction.

clear all
close all
% warning off all

%% IMPORT INPUT DATA
clear all


namefileinp = 'Input_v9.xlsx';
sheetinp = 'Inp_gruop_Rotter';%
sheetinp2 = 'Inp_plastic_Rotter';%
% switch_soil=1;
switch_soil=9000;
    rho_deep=-1;
    switch_shape=1;
switch_limf=1; %1 alpha, 2 beta, 3 CE 4Tabular
switch_displfound=2; %1 displacment 2 nondisplacement
switch_ep=1;   %1 ep or nep, 0 el
switch_nep=0;  %1 ep, 0 nep
% % % switch_uip=0;%to remove sliders in the nonlinear analysis
switch_headspring=0;
switch_cap=0; % 0 none, 1 kinematic constraint, 2 beam grid
coeffux=0; % 1-eps=1, 0-no. GF soil horizontal movements
coeffuz=0.00001; % 1-eps=1, 0-no. GF soil vertical movements
%=========================================
% switch_chenflex=600;
switch_type=0; % 0 uniform 1 multylayer
TPNum_global=3;  % 1 point load, 2 uniform circular, 3 ring load 4 uniform load for the vertcal settlement of base due to force on base; ring load is used for other case
%============================================
switch_sliderx=0; %activate the plasticity in the horizontal direction

switch_headdeg=1;

% if switch_limf=4 introduce tabular limit pressures
% xstresslim_vector_tabular=ones(11,1);
% xstresslim_vector_tabular(1:6)=6*6*1000;
% xstresslim_vector_tabular(7:end)=500*1000;
% zstresslim_vector_tabular=500*1000*ones(11,1);
% 
% zstresslim_vector_tabular=[
% 4800
% 4800
% 4800
% 4800
% 4800
% 4800
% 5000000
% 5000000
% 5000000
% 5000000
% 5000000];
% 
% 
% 
% xstresslim_vector_tabular=[
% 2000
% 2000
% 2000
% 27000
% 27000
% 27000
% 5000000
% 5000000
% 5000000
% 5000000
% 5000000];

switch_negative_fric=1; % half of the greenfield displacement history is to apply a linear settlement of surface value settlement_history_surface
settlement_history_surface=0.2;

switch_modfc=0; %calculation of modification factors

run('Analisi_staticinput_cont_v9')
run('Analisi_static_v13MSCA_Clash')
script_dz
%%

figure(11);hold on;
scale=10;
x_global_fin=x_global+Uxpg_fr_i_ep_deltaT*scale;
z_global_fin=z_global+Uzpg_fr_i_ep_deltaT*scale;
plot(x_global_fin,z_global_fin,'bo')
plot(x_global,z_global,'k+')

x_global_fin=x_global+Uxpg_fr_i_el_deltaT*scale;
z_global_fin=z_global+Uzpg_fr_i_el_deltaT*scale;
plot(x_global_fin,z_global_fin,'rx')
plot(x_global,z_global,'k+')

set ( gca, 'ydir', 'reverse' )
xlabel('X (m)')
ylabel('Z (m)')
pbaspect([1 1 1])
xlim([-15 15])
ylim([0 30])
%% Chen's example
% close all
[Ned,M_rd] = CosenzaCicular(d,0.02);

plot(Ned,M_rd,'-');hold on;
plot(Ned,-M_rd,'-')

plot(Npg_fr_i_el_T,Mpg_fr_i_el_T,'x'); hold on;

figure(2)

% plot(Mpg_fr_i_ep_T,z_global,'k--'); hold on;

% plot(Mpg_fr_i_ep_P,z_global,'r--'); hold on;

plot(Npg_fr_i_ep_P,z_global,'r-'); hold on;
plot(Npg_fr_i_el_P,z_global,'r:o'); hold on;
plot(Npg_fr_i_ep_T,z_global,'k-'); hold on;



set(gca, 'YDir','reverse')
figure(3)
subplot(1,2,1)
% plot(stresszpg_fr_i_ep_P(1:end-1),z_global(1:end-1),'r-'); hold on;
plot(stresszpg_fr_i_ep_T(1:end-1),z_global(1:end-1),'r--'); hold on;

plot(stresszpg_fr_i_el_T(1:end-1),z_global(1:end-1),'rx'); hold on;

% plot(stresszpg_fr_i_el_P(1:end-1)./flim1_vector(1:end-1),z_global(1:end-1),'r-x'); hold on;
% plot(stresszpg_fr_i_el_T(1:end-1)./flim1_vector(1:end-1),z_global(1:end-1),'r--x'); hold on;
set(gca, 'YDir','reverse')
subplot(1,2,2)
% plot(stressxpg_fr_i_ep_P,z_global,'k-'); hold on;
plot(stressxpg_fr_i_ep_T,z_global,'k--'); hold on;

% plot(stressxpg_fr_i_el_P./cu/9,z_global,'k-x'); hold on;
plot(stressxpg_fr_i_el_T,z_global,'kx'); hold on;

set(gca, 'YDir','reverse')





Cap_result=[x_centre/ht,U_cap_i_master,U_cap_h_master];
figure(199)
subplot(1,2,1)
plot(Sz(:,1),z,'gs-');set(gca,'Ydir','reverse');hold on;
plot(Uzpg_free(:,1),z,'r-');set(gca,'Ydir','reverse');hold on;
% plot(Uzpg_free(:,2),z,'k-');
plot(Uzpg_cap_i(:,1),z,'rx');
% plot(Uzpg_cap_i(:,2),z,'kx-');
plot(Uzpg_fr_i_el_deltaT(:,1),z,'ro');
plot(Uzpg_fr_i_ep_deltaT(:,1),z,'rs');

subplot(1,2,2)
plot(Sx(:,1),z,'gs--');set(gca,'Ydir','reverse');hold on;
plot(Uxpg_free(:,1),z,'k-');set(gca,'Ydir','reverse');hold on;
% plot(Uxpg_free(:,2),z,'k-');
plot(Uxpg_cap_i(:,1),z,'kx');
% plot(Uxpg_cap_i(:,2),z,'kx-');
plot(Uxpg_fr_i_ep_deltaT(:,1),z,'ko');
plot(Uxpg_fr_i_el_deltaT(:,1),z,'ks');

%%
% close all
figure(1099); hold on;
subplot(1,2,1)
plot(Xi_nod_base,Sz(1,:),'gs-');
set(gca,'Ydir','reverse');hold on;
% plot(Xi_nod_base,Uzpg_free_el_deltaT(1,:),'r-');
plot(Xi_nod_base,Uzpg_free_ep_deltaT(1,:),'r--');
% plot(Uzpg_free(:,2),z,'k-');
% plot(Xi_nod_base,Uzpg_cap_i(1,:),'rx');
% plot(Uzpg_cap_i(:,2),z,'kx-');
% plot(Xi_nod_base,Uzpg_fr_i_el_deltaT(1,:),'ro');
plot(Xi_nod_base,Uzpg_fr_i_ep_deltaT(1,:),'rs');
legend('Greenfield - surface','Free Head - Elastic','Free Head - Nonlinear','Active 2st Building - Elastic','Active 2st Building - Nonlinear','location', 'northoutside')
legend('Greenfield - surface','Free Head - Nonlinear','Active 2st Building - Nonlinear','location', 'northoutside')

ylabel('Settlement (m)')
xlabel('Offset X (m)')

subplot(1,2,2); hold on;
scale=300;
x_global_fin=x_global+Uxpg_fr_i_ep_deltaT*scale;
z_global_fin=z_global+Uzpg_fr_i_ep_deltaT*scale;
plot(x_global,z_global,'k+')
plot(x_global_fin,z_global_fin,'bo')


x_global_fin=x_global+Uxpg_free_ep_deltaT*scale;
z_global_fin=z_global+Uzpg_free_ep_deltaT*scale;
plot(x_global_fin,z_global_fin,'rx')


set ( gca, 'ydir', 'reverse' )
xlabel('X (m)')
ylabel('Z (m)')
pbaspect([1 1 1])
xlim([-15 15])
ylim([0 30])

title('Foundation deformed shape: Nonlinear analyses')%,'Foundation deformed shape','location', 'northoutside'

%%
figure(2122)
plot(degscalar/degscalar(1)*100)
xlabel('time')
ylabel('head stiffness %')
legend('reduction stiffness selected piles','location', 'northoutside')
%%
% plot(Npg_free(:,1)/1000,z,'r-');set(gca,'Ydir','reverse');set(gca,'Xdir','reverse');hold on;
% plot(Npg_cap_i(:,1)/1000,z,'rx');set(gca,'Ydir','reverse');set(gca,'Xdir','reverse');hold on;

% plot(Npg_free(:,2)/1000,z,'k-');set(gca,'Ydir','reverse');set(gca,'Xdir','reverse');hold on;
% plot(Npg_cap_i(:,2)/1000,z,'kx-');set(gca,'Ydir','reverse');set(gca,'Xdir','reverse');hold on;


plot(Mpg_free(:,1)/1000,z,'r-');set(gca,'Ydir','reverse');set(gca,'Xdir','reverse');hold on;
plot(Mpg_cap_i(:,1)/1000,z,'rx');set(gca,'Ydir','reverse');set(gca,'Xdir','reverse');hold on;

% plot(Mpg_free(:,2)/1000,z,'k-');set(gca,'Ydir','reverse');set(gca,'Xdir','reverse');hold on;
% plot(Mpg_cap_i(:,2)/1000,z,'kx-');set(gca,'Ydir','reverse');set(gca,'Xdir','reverse');hold on;


%%
close all
plot(stresspg_fr_i_el_deltaT(1:end-1),z_global(1:end-1),'k');hold on
plot(stresspg_fr_i_ep_deltaT(1:end-1),z_global(1:end-1),'r');hold on
%%
close all


figure(1)
% plot(F_soil_fr_i_el_P_z(:,end),'b-');hold on
plot(F_soil_fr_i_ep_P_z(:,end),'b--');hold on
% plot(F_soil_fr_i_el_P_x(:,end),'b-');hold on
% plot(F_soil_fr_i_ep_P_x(:,end),'b--');hold on

% plot(F_soil_fr_i_el_T_z(:,end),'k-');hold on
plot(F_soil_fr_i_ep_T_z(:,end),'k--');hold on
% plot(F_soil_fr_i_el_T_x(:,end),'r-');hold on
% plot(F_soil_fr_i_ep_T_x(:,end),'r--');hold on

figure(2)
plot(vect_Utotf_fr_i_el_P_z(:,end),'b-');hold on
plot(vect_Utotf_fr_i_ep_P_z(:,end),'b--')
% plot(vect_Utotf_fr_i_el_P_x(:,end),'b-')
% plot(vect_Utotf_fr_i_ep_P_x(:,end),'b--')

plot(vect_Utotp_fr_i_el_P_z(:,end),'bo');hold on
plot(vect_Utotp_fr_i_ep_P_z(:,end),'bo')
% plot(vect_Utotp_fr_i_el_P_x(:,end),'bo')
% plot(vect_Utotp_fr_i_ep_P_x(:,end),'bo')

plot(vect_Utotf_fr_i_el_deltaT_z(:,end),'g-');hold on
plot(vect_Utotf_fr_i_ep_deltaT_z(:,end),'g--')
% plot(vect_Utotf_fr_i_el_deltaT_x(:,end),'r-')
% plot(vect_Utotf_fr_i_ep_deltaT_x(:,end),'r--')

plot(vect_Utotf_fr_i_el_T_z(:,end),'k-');hold on
plot(vect_Utotf_fr_i_ep_T_z(:,end),'k--')
% plot(vect_Utotf_fr_i_el_deltaT_x(:,end),'r-')
% plot(vect_Utotf_fr_i_ep_deltaT_x(:,end),'r--')

% plot(vect_Utotp_fr_i_el_T_z(:,end),'ko');hold on
plot(vect_Utotp_fr_i_ep_T_z(:,end),'ko')
% plot(vect_Utotp_fr_i_el_T_x(:,end),'ro')
% plot(vect_Utotp_fr_i_ep_T_x(:,end),'ro')





figure(3)
h(11)=plot(Xi_nod_base,vect_Utotf_fr_i_el_P_z(1:nnodes:end,end),'ob-');hold on
h(12)=plot(Xi_nod_base,vect_Utotf_fr_i_ep_P_z(1:nnodes:end,end),'ob--');

h(21)=plot(Xi_nod_base,vect_Utotf_fr_i_el_deltaT_z(1:nnodes:end,end),'og-');hold on
h(22)=plot(Xi_nod_base,vect_Utotf_fr_i_ep_deltaT_z(1:nnodes:end,end),'og--');

h(31)=plot(Xi_nod_base,vect_Utotf_fr_i_el_T_z(1:nnodes:end,end),'ok-');hold on
h(32)=plot(Xi_nod_base,vect_Utotf_fr_i_ep_T_z(1:nnodes:end,end),'ok--');

title(['Piles Vlt=',num2str(Vltp,2),', SFo=',num2str(SF0,2),''],'FontSize',9)

xlabel('Offset (m)','FontSize',9);% 
ylabel('Settlement (m)','FontSize',9);

Legend1 =legend( [h(11),h(12),h(21),h(22),h(31),h(32)],...
    {'Loads - EL','Loads - EP','Tunn - EL','Tunn - EP','Tunn + Loads - EL','Tunn + Loads - EP'},...
    'Location','eastoutside','Fontsize',9,'EdgeColor',[1 1 1]);%,'PlotBoxAspectRatio',ratio_leg21'Position',pos_leg_1','PlotBoxAspectRatioMode','manual'
legend('boxon')
set(gca,'Ydir','reverse')
% title(%


figure(4)
h(101)=plot(Npg_fr_i_ep(2:end,1),z(2:end),'b-');hold on

title(['Piles Vlt=',num2str(Vltp,2),', SFo=',num2str(SF0,2),''],'FontSize',9)

xlabel('Axial force (N)','FontSize',9);% 
ylabel('Depth (m)','FontSize',9);

% Legend1 =legend( [h(11),h(12),h(21),h(22),h(31),h(32)],...
%     {'Loads - EL','Loads - EP','Tunn - EL','Tunn - EP','Tunn + Loads - EL','Tunn + Loads - EP'},...
%     'Location','eastoutside','Fontsize',9,'EdgeColor',[1 1 1]);%,'PlotBoxAspectRatio',ratio_leg21'Position',pos_leg_1','PlotBoxAspectRatioMode','manual'
% legend('boxon')
set(gca,'Ydir','reverse')
% title(%



figure(5)
h(11)=plot(0,vect_Utotf_free_ep_P_z(1,end),'ob-');hold on

h(21)=plot(Vltp_vect(2:end),vect_Utotf_free_ep_deltaT_z(1,:),'og-');hold on
h(21)=plot([0,Vltp_vect(end)],[0,Uzpg_free(1,1)],'g:');hold on

h(31)=plot(Vltp_vect(2:end),vect_Utotf_free_ep_T_z(1,:),'ok-');hold on
axis([0 5 0 0.2])
title(['Piles Phead=',num2str(max(P_el)/1000,2),'kN'],'FontSize',9)

xlabel('Vlt (%)','FontSize',9);% 
ylabel('Settlement (m)','FontSize',9);

Legend1 =legend( [h(11),h(12),h(21),h(22),h(31),h(32)],...
    {'Loads - EL','Loads - EP','Tunn - EL','Tunn - EP','Tunn + Loads - EL','Tunn + Loads - EP'},...
    'Location','eastoutside','Fontsize',9,'EdgeColor',[1 1 1]);%,'PlotBoxAspectRatio',ratio_leg21'Position',pos_leg_1','PlotBoxAspectRatioMode','manual'
legend('boxon')
set(gca,'Ydir','reverse')
% title(%

stop

%% Deep excavation soil999
close all

figure(6)

h(1)=plot(vect_Utotf_fr_i_ep_deltaT_z(1:end,end),z,'or-');hold on
h(2)=plot(vect_Utotf_fr_i_el_deltaT_z(1:end,end),z,'ok-');
h(3)=plot(vect_Ugf_z(1:end,end),z,'og-');

h(1)=plot(vect_Utotf_fr_i_ep_deltaT_z(1:end,end),z,'or-');hold on
h(2)=plot(vect_Utotf_fr_i_el_deltaT_z(1:end,end),z,'ok-');
h(3)=plot(vect_Ugf_z(1:end,end),z,'og-');

% axis([0 5*So 0 max(z)])
title(['SF=',num2str(max(SF0),2),''],'FontSize',9)

ylabel('Vlt (%)','FontSize',9);% 
xlabel('Depth (m)','FontSize',9);
set(gca,'Ydir','reverse')
Legend1 =legend( [h(1),h(2),h(3)],...
    {'EP','EL','GF'},...
    'Location','eastoutside','Fontsize',9,'EdgeColor',[1 1 1]);%,'PlotBoxAspectRatio',ratio_leg21'Position',pos_leg_1','PlotBoxAspectRatioMode','manual'
legend('boxon')

% title(%

Dz=[0.0114385965336643];
bb=[-deltaS/Dz,1/SF0,-(vect_Utotf_fr_i_ep_deltaT_z(1,end)-So)/So];
aa=[P_el(3),vect_Utotf_fr_i_el_P_z(1,end),vect_Utotf_free_ep_P_z(1,end)];

%% check hyperbolic stiffness degradation

 
close all
%loading
P_elvect=P_el(3)/(num_incr_P):P_el(3)/(num_incr_P):P_el(3);
%loading and unloaidng
% P_elvect1=[P_el(3)/(num_incr_P):P_el(3)/(num_incr_P):P_el(3)/2];
% P_elvect2=[P_el(3)/2:-P_el(3)/(num_incr_P):0];
% P_elvect=[P_elvect1,P_elvect2];


Pu=P_el(3)*SF0;
Kel=P_elvect(2)/vect_Utotf_free_ep_P(3,2);
Dz=Pu/Kel;


plot(P_elvect/10^6,vect_Utotf_free_ep_P(3,:)/d*100,'b');hold on;
% plot(P_elvect1,vect_Utotf_free_ep_P(3,1:length(P_elvect1))/d*100,'r');hold on;
plot([0 Kel*vect_Utotf_free_ep_P(3,end)]/10^6,[0 vect_Utotf_free_ep_P(3,end)]/d*100,'-om')
plot([0 Kel*vect_Utotf_free_ep_P(3,end)]/10^6,[Dz Dz]/d*100,'-k')

%%

figure(199)

plot(vect_Utotf_fr_i_el_P_z(1,:));

%%
%%
close all
figure(111)
subplot(2,3,1)
vv=Uffz(:,:,end);
plot(X,Uzpg_fr_i_el_deltaT(1,:),'k-x');hold on;
plot(X,Uzpg_fr_i_ep_deltaT(1,:),'k:o');hold on;
plot(X,Uzpg_free_el_deltaT(1,:),'r-x');hold on;
plot(X,Uzpg_free_ep_deltaT(1,:),'r:o');hold on;
plot(X,Sz(1,:),'b:o');
plot(X,Sz(end,:),'b-o');
title('Pile heads')
xlabel('offset from wall (m)');ylabel('excavation induced settlement (m)')
legend('pile head with beam - EL','pile head with beam - EP','pile head freehead - EL','pile head freehead - EP', 'GF - head', 'GF - tip', 'location','best')
set(gca, 'YDir','reverse')

subplot(2,3,2)
plot(X,Uxpg_fr_i_el_deltaT(1,:),'k-x');hold on;
plot(X,Uxpg_fr_i_ep_deltaT(1,:),'k:o');hold on;
plot(X,Uxpg_free_el_deltaT(1,:),'r-x');hold on;
plot(X,Uxpg_free_ep_deltaT(1,:),'r:o');hold on;
plot(X,Sx(1,:),'b:o');
plot(X,Sx(end,:),'b-o');
title('Pile heads')
xlabel('offset from wall (m)');ylabel('excavation induced horizontal movement (m)')
legend('pile head with beam - EL','pile head with beam - EP','pile head freehead - EL','pile head freehead - EP', 'GF - head', 'GF - tip', 'location','best')
set(gca, 'YDir','reverse')

subplot(2,3,3)
plot(Uzpg_fr_i_el_deltaT(:,1),z,'k-x');hold on;
plot(Uzpg_fr_i_ep_deltaT(:,1),z,'k:o');hold on;
plot(Uzpg_free_el_deltaT(:,1),z,'r-x');hold on;
plot(Uzpg_free_ep_deltaT(:,1),z,'r:o');hold on;
plot(Sz(:,1),z,'b-');
title('Uz - first pile')
ylabel('depth (m)');xlabel('Settlement')
legend('pile with beam - EL','pile with beam - EP','pile freehead - EP','pile freehead - EP','GF', 'location','best')
set(gca, 'YDir','reverse')

subplot(2,3,4)
plot(Uzpg_fr_i_el_deltaT(:,2),z,'k-x');hold on;
plot(Uzpg_fr_i_ep_deltaT(:,2),z,'k:o');hold on;
plot(Uzpg_free_el_deltaT(:,2),z,'r-x');hold on;
plot(Uzpg_free_ep_deltaT(:,2),z,'r:o');hold on;
plot(Sz(:,2),z,'b-');
title('Uz - second pile')
ylabel('depth (m)');xlabel('Settlement')
legend('pile with beam - EL','pile with beam - EP','pile freehead - EP','pile freehead - EP','GF', 'location','best')
set(gca, 'YDir','reverse')

subplot(2,3,5)
plot(Uxpg_fr_i_el_deltaT(:,1),z,'k-x');hold on;
plot(Uxpg_fr_i_ep_deltaT(:,1),z,'k:o');hold on;
plot(Uxpg_free_el_deltaT(:,1),z,'r-x');hold on;
plot(Uxpg_free_ep_deltaT(:,1),z,'r:o');hold on;
plot(Sx(:,1),z,'b-');
title('Ux - first pile')
ylabel('depth (m)');xlabel('hor movement')
legend('pile with beam - EL','pile with beam - EP','pile freehead - EP','pile freehead - EP','GF', 'location','best')
set(gca, 'YDir','reverse')

subplot(2,3,6)
plot(Uxpg_fr_i_el_deltaT(:,2),z,'k-x');hold on;
plot(Uxpg_fr_i_ep_deltaT(:,2),z,'k:o');hold on;
plot(Uxpg_free_el_deltaT(:,2),z,'r-x');hold on;
plot(Uxpg_free_ep_deltaT(:,2),z,'r:o');hold on;
plot(Sx(:,2),z,'b-');
title('Ux - second pile')
ylabel('depth (m)');xlabel('hor movement')
legend('pile with beam - EL','pile with beam - EP','pile freehead - EP','pile freehead - EP','GF', 'location','best')
set(gca, 'YDir','reverse')

%%%%%%%%%%%%%%%%%%%%%%
figure(200)
subplot(2,3,1)
vv=Uffz(:,:,end);
plot(X,Npg_fr_i_ep_deltaT(1,:),'k-x');hold on;
plot(X,Npg_fr_i_ep_deltaT(1,:),'k:o');hold on;
plot(X,Npg_free_ep_deltaT(1,:),'r-x');hold on;
plot(X,Npg_free_ep_deltaT(1,:),'r:o');hold on;
title('Pile heads')
xlabel('offset from wall (m)');ylabel('excavation induced axial force (Nm)')
legend('pile with beam - EL','pile with beam - EP','pile freehead - EP','pile freehead - EP', 'location','best')
set(gca, 'YDir','reverse')

subplot(2,3,2)
plot(X,Mpg_fr_i_ep_deltaT(1,:),'k-x');hold on;
plot(X,Mpg_fr_i_ep_deltaT(1,:),'k:o');hold on;
plot(X,Mpg_free_ep_deltaT(1,:),'r-x');hold on;
plot(X,Mpg_free_ep_deltaT(1,:),'r:o');hold on;
title('Pile heads')
xlabel('offset from wall (m)');ylabel('excavation induced bending moment (Nm)')
legend('pile with beam - EL','pile with beam - EP','pile freehead - EP','pile freehead - EP', 'location','best')
set(gca, 'YDir','reverse')

subplot(2,3,3)
plot(Npg_fr_i_ep_deltaT(:,1),z,'k-x');hold on;
plot(Npg_fr_i_ep_deltaT(:,1),z,'k:o');hold on;
plot(Npg_free_ep_deltaT(:,1),z,'r-x');hold on;
plot(Npg_free_ep_deltaT(:,1),z,'r:o');hold on;
title('N - first pile')
ylabel('depth (m)');xlabel('Axial force')
legend('pile with beam - EL','pile with beam - EP','pile freehead - EP','pile freehead - EP', 'location','best')
set(gca, 'YDir','reverse')

subplot(2,3,4)
plot(Npg_fr_i_ep_deltaT(:,2),z,'k-x');hold on;
plot(Npg_fr_i_ep_deltaT(:,2),z,'k:o');hold on;
plot(Npg_free_ep_deltaT(:,2),z,'r-x');hold on;
plot(Npg_free_ep_deltaT(:,2),z,'r:o');hold on;
title('N - second pile')
ylabel('depth (m)');xlabel('Axial force')
legend('pile with beam - EL','pile with beam - EP','pile freehead - EP','pile freehead - EP', 'location','best')
set(gca, 'YDir','reverse')

subplot(2,3,5)
plot(Mpg_fr_i_ep_deltaT(:,1),z,'k-x');hold on;
plot(Mpg_fr_i_ep_deltaT(:,1),z,'k:o');hold on;
plot(Mpg_free_ep_deltaT(:,1),z,'r-x');hold on;
plot(Mpg_free_ep_deltaT(:,1),z,'r:o');hold on;
title('M - first pile')
ylabel('depth (m)');xlabel('Bending moment')
legend('pile with beam - EL','pile with beam - EP','pile freehead - EP','pile freehead - EP', 'location','best')
set(gca, 'YDir','reverse')

subplot(2,3,6)
plot(Mpg_fr_i_ep_deltaT(:,2),z,'k-x');hold on;
plot(Mpg_fr_i_ep_deltaT(:,2),z,'k:o');hold on;
plot(Mpg_free_ep_deltaT(:,2),z,'r-x');hold on;
plot(Mpg_free_ep_deltaT(:,2),z,'r:o');hold on;
title('M - second pile')
ylabel('depth (m)');xlabel('Bending moment')
legend('pile with beam - EL','pile with beam - EP','pile freehead - EP','pile freehead - EP', 'location','best')
set(gca, 'YDir','reverse')

%%
% 
% stop
% %% PLOT
% 
% 
% figure(1)
% %set(1, 'Units','centimeters', 'Position',[0 0 18 6])
% 
% 
% hold on
% plot(X(1,:),Sz(1,:)*1000,':k','LineWidth',3)
% plot(X(1,:),Uzpg_free(1,:)*1000,'r','LineWidth',1.5)
% plot(X(1,:),Uzpg_avg(1,:)*1000,'k','LineWidth',1.5) %WITHOUT BASE
% plot(X(1,:),Uzpg_rig_base(1,:)*1000,'g','LineWidth',2.5)%WITH BASE
% plot(X(1,:),Uzpg_fr_i(1,:)*1000,'b','LineWidth',1.5)
% plot(X(1,:),Uzpg_fr_h(1,:)*1000,'m','LineWidth',1.5)
% plot(X(1,:),Uzpg_fr_vspring(1:npali)*1000,':y','LineWidth',3)
% %title(['Lateral Deflection (mm), Front Pile'],'FontSize',11.0,'FontName','Cambria')
% ylabel('Settlement (mm)','FontSize',11.0,'FontName','Cambria')
% xlabel('z respect tunnel CL (m)','FontSize',11.0,'FontName','Cambria')
% set(gca,'YDir','reverse');
% %         ylim([-5 50])
% %axis([-9 0 0 25])
% legend({'Greenfield','Tunnel-pile group','Tunnel-pile avg (no base)','Tunnel-pile avg (with base)','Tunnel-pile-frame (fixed pile head)','Tunnel-pile-frame (hinged pile head)','Tunnel-spring-frame (hinged pile head)'},'Location','Northwest','FontSize',8,'FontName','Cambria')
% legend('boxoff')
% 
% 
% figure(2)
% %set(1, 'Units','centimeters', 'Position',[0 0 18 6])
% 
% 
% hold on
% plot(X(1,:),Uzpg_fr_i(1,:)*1000,'b','LineWidth',1.5)
% plot(X(1,:),Uzpg_fr_h(1,:)*1000,'m','LineWidth',1.5)
% plot(X(1,:),Uzpg_fr_vspring(1:npali)*1000,':y','LineWidth',3)
% %title(['Lateral Deflection (mm), Front Pile'],'FontSize',11.0,'FontName','Cambria')
% ylabel('Settlement (mm)','FontSize',11.0,'FontName','Cambria')
% xlabel('z respect tunnel CL (m)','FontSize',11.0,'FontName','Cambria')
% set(gca,'YDir','reverse');
% %         ylim([-5 50])
% %axis([-9 0 0 25])
% legend({'Tunnel-pile-frame (fixed pile head)','Tunnel-pile-frame (hinged pile head)','Tunnel-spring-frame (hinged pile head)'},'Location','Northwest','FontSize',8,'FontName','Cambria')
% legend('boxoff')
% 
% 
% figure(3)
% %set(1, 'Units','centimeters', 'Position',[0 0 18 6])
% 
% hold on
% plot(X(1,:),Uzpg_free(1,:)*1000,'r','LineWidth',1.5)
% plot(X(1,:),Uzpg_avg(1,:)*1000,'k','LineWidth',1.5) %WITHOUT BASE
% plot(X(1,:),Uzpg_rig_base(1,:)*1000,'g','LineWidth',2.5)%WITH BASE
% plot(X(1,:),Sz(1,:)*1000,':k','LineWidth',3)
% 
% %title(['Lateral Deflection (mm), Front Pile'],'FontSize',11.0,'FontName','Cambria')
% ylabel('Settlement (mm)','FontSize',11.0,'FontName','Cambria')
% xlabel('z respect tunnel CL (m)','FontSize',11.0,'FontName','Cambria')
% set(gca,'YDir','reverse');
% %         ylim([-5 50])
% %axis([-9 0 0 25])
% legend({'Tunnel-pile group','Tunnel-pile avg (no base)','Tunnel-pile avg (with base)','Greenfield'},'Location','Northwest','FontSize',8,'FontName','Cambria')
% legend('boxoff')
% 
% figure(4)
% %set(1, 'Units','centimeters', 'Position',[0 0 18 6])
% 
% 
% hold on
% plot(X(1,:),Uzpg_free(1,:)*1000,'r','LineWidth',1.5)
% plot(X(1,:),Uzpg_fr_i(1,:)*1000,'b','LineWidth',1.5)
% plot(X(1,:),Uzpg_fr_vspring(1:npali)*1000,':y','LineWidth',3)
% % % % %         plot(X(1,:),Uzpg_fr_i_c_eq_beam(1,:)*1000,':b','LineWidth',3)
% %title(['Lateral Deflection (mm), Front Pile'],'FontSize',11.0,'FontName','Cambria')
% ylabel('Settlement (mm)','FontSize',11.0,'FontName','Cambria')
% xlabel('z respect tunnel CL (m)','FontSize',11.0,'FontName','Cambria')
% set(gca,'YDir','reverse');
% %         ylim([-5 50])
% %axis([-9 0 0 25])
% legend({'Tunnel-pile group','Tunnel-pile-frame (fixed pile head)','Tunnel-spring-frame (hinged pile head)','Tunnel-pile-frame with eq beam (fixed pile head)'},'Location','Northwest','FontSize',8,'FontName','Cambria')
% legend('boxoff')
% 
% hold off
% 
% Excelz=[Sz(1,:)',Uzpg_free(1,:)',Uzpg_fr_i(1,:)',Uzpg_fr_h(1,:)',Uzpg_fr_vspring']*1000;
% Excelx=[Sx(1,:)',Uxpg_free(1,:)',Uxpg_fr_i(1,:)',Uxpg_fr_h(1,:)',Uxpg_fr_h(1,:)'*0]*1000;
% Excelr2=[0*Sz(1,:)',Ur2pg_free(1,:)',Ur2pg_fr_i(1,:)',Ur2pg_fr_h(1,:)',Uxpg_fr_h(1,:)'*0];
% 
% %% PLOT
% up_bond_x=[1e-6,1e-3,1,1e2]*10; %*10 to correct EI/m
% up_bond_y=[1.2,1.2,0,0];
% 
% dw_bond_x=[1e-6,1e-4,5e-02,1e2]*10;%*10 to correct EI/m
% dw_bond_y=[1,1,0,0];
% 
% 
% size_picture=[ 2 2 10 10];
% size_ticklength=2;
% mrkrsize=3.0;
% figure(110); clf;
% set(gcf, 'Units','centimeters','PaperUnits','centimeters','PaperPositionMode','auto','Position',size_picture)
% set(gca,'FontName','Times New Roman','FontSize',9,'linewidth',1.0);
% set(gca,'ticklength',size_ticklength*get(gca,'ticklength'))
% 
% % h(1)=semilogx(up_bond_x,up_bond_y,'b-');
% % hold on
% % h(2)=semilogx(dw_bond_x,dw_bond_y,'b-');
% % h(3)=semilogx(data_red_factor(:,18),data_red_factor(:,23),'ko');%\rho_{sag} vs M^{DR,sag}-Fixed Head
% % h(4)=semilogx(data_red_factor(:,19),data_red_factor(:,24),'k^');%\rho_{hog} vs M^{DR,hog}-Fixed Head
% % h(5)=line([data_red_factor(:,18),data_red_factor(:,18)]',[data_red_factor(:,21),data_red_factor(:,23)]','color','k');%red_factor_reldefl_sag_free
% % h(5)=line([data_red_factor(:,19),data_red_factor(:,19)]',[data_red_factor(:,22),data_red_factor(:,24)]','color','k');%red_factor_reldefl_hog_free
% % h(7)=semilogx([1e-5 1e5],[1 1],'g-');%red_factor_reldefl_hog_free
% % title(['Modification Factors vs Relative Building Stiffness'],'FontSize',10)
% % xlabel('\rho_{sag};\rho_{hog} (-)','FontSize',10);
% % ylabel('M^{DR,sag};M^{DR,hog} (-)','FontSize',10);
% % legend (h([7,1,3,4]),{'Greenfield','Mair(2013)','SAG-Fixed Head','HOG-Fixed Head'},'Location','Northeast')
% % 


