%% check hyperbolic stiffness degradation
 figure(1000);hold on;
% close all
%loading
P_elvect=P_el(3)/(num_incr_P):P_el(3)/(num_incr_P):P_el(3);

P_elvect=P_el_vect(3,:);


%loading and unloaidng
% P_elvect1=[P_el(3)/(num_incr_P):P_el(3)/(num_incr_P):P_el(3)/2];
% P_elvect2=[P_el(3)/2:-P_el(3)/(num_incr_P):0];
% P_elvect=[P_elvect1,P_elvect2];


Pu=P_el(3)*SF0;
Kel=P_elvect(1)/vect_Utotf_free_ep_P(3,1);
Dz=Pu/Kel;


% plot(P_elvect/10^6,vect_Utotf_free_ep_P(3,:)/d*100,'b');hold on;
% % plot(P_elvect1,vect_Utotf_free_ep_P(3,1:length(P_elvect1))/d*100,'r');hold on;
% plot([0 Kel*vect_Utotf_free_ep_P(3,end)]/10^6,[0 vect_Utotf_free_ep_P(3,end)]/d*100,'-om')
% plot([0 Kel*vect_Utotf_free_ep_P(3,end)]/10^6,[Dz Dz]/d*100,'-k')
% figure(2)
plot(P_elvect,vect_Utotf_free_ep_P(3,:),'-b');hold on;
plot(P_elvect,vect_Utotf_fr_i_ep_P(3,:),'--b');hold on;
xlabel('Force (N)')
ylabel('uz (m)')
% title(['e_x/d= ' num2str((ecc_Ncap_x/d)) '  SF0=' num2str((SF0)) ''])
plot([0 P_elvect(end)],[0,vect_Utotf_fr_i_el_P(3,end)],'sc','markersize',10);hold on;
% plot(P_elvect,vect_Utotf_free_el_P(3,:),'-k');hold on;



% plot(P_elvect1,vect_Utotf_free_ep_P(3,1:length(P_elvect1))/d*100,'r');hold on;
plot([0 Kel*vect_Utotf_free_ep_P(3,end)],[0 vect_Utotf_free_ep_P(3,end)],'-om')
% plot([0 Kel*vect_Utotf_free_ep_P(3,end)],[Dz Dz],'-k')


figure(2000)
plot(P_elvect,vect_Utotf_free_ep_P(5,:),'-b');hold on;
plot(P_elvect,vect_Utotf_fr_i_ep_P(5,:),'--b');hold on;
xlabel('Force (N)')
ylabel('Rotation_y (rad)')
% title(['e_x/d= ' num2str((ecc_Ncap_x/d)) '  SF0=' num2str((SF0)) ''])

figure(3000)
plot(P_elvect,vect_Utotf_free_ep_P(1,:),'-b');hold on;
plot(P_elvect,vect_Utotf_fr_i_ep_P(1,:),'--b');hold on;
xlabel('Force (N)')
ylabel('ux (m)')
% title(['e_x/d= ' num2str((ecc_Ncap_x/d)) '    SF0=' num2str((SF0)) ''])


%%
figure(10);hold on;
scale=200;
x_global_fin=x_global+Uxpg_fr_i_ep_P*scale;
z_global_fin=z_global+Uzpg_fr_i_ep_P*scale;
plot(x_global_fin,z_global_fin,'bo')
plot(x_global,z_global,'k+')

x_global_fin=x_global+Uxpg_fr_i_el_P*scale;
z_global_fin=z_global+Uzpg_fr_i_el_P*scale;
plot(x_global_fin,z_global_fin,'rx')
plot(x_global,z_global,'k+')

set ( gca, 'ydir', 'reverse' )
xlabel('X (m)')
ylabel('Z (m)')
pbaspect([1 1 1])
xlim([-15 15])
ylim([0 30])
% title(['e_x/d= ' num2str((ecc_Ncap_x/d)) '  scale factor=' num2str((scale)) '  SF0=' num2str((SF0)) ''])

%%
% Table_Brian_EP=[npali*P_elvect/1000;vect_Utotf_fr_i_ep_P(3,:);vect_Utotf_fr_i_ep_P(1,:);vect_Utotf_fr_i_ep_P(5,:);vect_Utotf_fr_i_ep_P(4,:)];
% Table_Brian_EP=-Table_Brian_EP';

% Table_Brian_EL=[P_el(3);P_el(5);ecc_Ncap_x/d;vect_Utotf_fr_i_el_P(3,end);vect_Utotf_fr_i_el_P(1,end);vect_Utotf_fr_i_el_P(5,end);vect_Utotf_fr_i_el_P(4,end)];
% Table_Brian_EL=-Table_Brian_EL';