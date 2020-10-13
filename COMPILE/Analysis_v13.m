%% Updates with respect to v12
% forig was added

%%


% clear all
% %
% 
% %==========================================================================
% %SWITCHERS
% %==========================================================================
% 
% switch_ep=0;    % 1 - yes   0 - no
% switch_soil=10; %1 clay, 2 sand
% 
% %%
% %==========================================================================
% %INPUT SCRIPT
% %==========================================================================
% Analysis_staticinput

%% ELABORATION INPUT DATA

%DISCRETIZZAZIONE DEL FOOTING

dx=diff(0:h_el_foot:B_X_foot);      %vettore della lunghezza degli elementi del singolo footing
nnodes_foot=B_X_foot./h_el_foot+1;
nfoot=length(X_foot_centr);         %numero dei footings
num_bay_foot=nfoot-1;
dimK=nnodes_foot*6;                 %dimensione matrice singolo footing
nelementi_foot=nnodes_foot-1;       %numero degli elementi del singolo palo


if nfoot==1
num_bay_raft=B_X_foot/spanx_frame;
integ=floor(num_bay_raft);
fract=num_bay_raft-integ;
if length(B_X_foot)>1
else
if fract>0
    stop
end
end
nnodes_bayraft=spanx_frame./h_el_foot+1;
nelementi_bayraft=nnodes_bayraft-1;
end




for i=1:nnodes_foot;
	for ip=1:nfoot;
       z_global_foot(i,ip)=Z_foot_centr(ip);
       x_global_foot(i,ip)=X_foot_centr(ip)-B_X_foot/2+(i-1)*h_el_foot;
       y_global_foot(i,ip)=Y_foot_centr(ip);
	end
end

for i=1:nfoot;
    znod_tot(1+(nelementi_foot+1)*(i-1):(nelementi_foot+1)+(nelementi_foot+1)*(i-1),1)=z_global_foot(1:nelementi_foot+1,i);
    xnod_tot(1+(nelementi_foot+1)*(i-1):1:(nelementi_foot+1)+(nelementi_foot+1)*(i-1),1)=x_global_foot(1:nelementi_foot+1,i) ;
    ynod_tot(1+(nelementi_foot+1)*(i-1):1:(nelementi_foot+1)+(nelementi_foot+1)*(i-1),1)=y_global_foot(1:nelementi_foot+1,i);
    Es_nod(1+(nelementi_foot+1)*(i-1):(nelementi_foot+1)+(nelementi_foot+1)*(i-1),1)=Es;
    ni_nod(1+(nelementi_foot+1)*(i-1):1:(nelementi_foot+1)+(nelementi_foot+1)*(i-1),1)=nis;
end
Gs_nod=Es_nod/2/(1+nis);
Zi_nod=znod_tot;
Xi_nod=xnod_tot;
Yi_nod=ynod_tot;

if range(diff(X_foot_centr)) ~= 0
h = msgbox('Transverse spacing of the footing needs to be uniform', 'Error','error');
stop
end

% % 
% % %% GREENFIELD DISPLACEMENT
% % % DESCRIPTIVE TEXT
% % 
% % if switch_soil==1
% %     [ Uffx,Uffz ]=u_LON(z_global_foot,x_global_foot,ht,Rt,epsilon,0,0.4999999);
% % else
% %     load out_u_reg_coeff_CD13_ID90.mat
% %     load out_u_reg_coeffx_CD13_ID90.mat
% %     id=0.9;
% %     [  Uffx,Uffz ]=u_SE(z_global_foot,x_global_foot,ht,Rt,epsilon,id,reg_ca,reg_cb,reg_c1,reg_c2,reg_c3,reg_c4,reg_c5,reg_c6,reg_cax,reg_cbx,reg_c1x,reg_c2x,reg_c3x,reg_c4x,reg_c5x,reg_c6x)  ;
% % end
% % Uffx=coeffux*Uffx;
% % Uffz=coeffuz*Uffz;
% % Uffy=zeros(nnodes_foot,nfoot);
% % 
% % Uffx(((z_global_foot-ht).^2+x_global_foot.^2).^0.5<=Rt)=0;
% % Uffz(((z_global_foot-ht).^2+x_global_foot.^2).^0.5<=Rt)=0;
% % 
% % %==========================================================================
% % %CARATTERISTICHE DELL'INPUT (SPOSTAMENTO DEL VINCOLO:"CEDIMENTO VINCOLARE")
% % %==========================================================================
% % %Quanto segue definisce l'input degli spostamenti nodali free-field. Il
% % %programma definisce la matrice Dfft che ha tante colonne per quanti sono
% % %gli istanti di analisi (comprese code di zeri) e tante righe quanti sono i
% % %gradi di libertà del singolo palo. I gradi di libertà sono ordinati a
% % %partire dal nodo in testa al palo secondo 
% % %  - spostamento orizzontale lungo x
% % %  - spostamento orizzontale lungo y
% % %  - spostamento orizzontale lungo z
% % %  - rotazione intorno a x
% % %  - rotazione intorno a y
% % %  - rotazione intorno a z
% % 
% %   
% % %Uff VECTOR AT NODES OF FINITE ELEMENTS
% % Dfft_gruppo=zeros(nfoot*nnodes_foot*6,1);
% % for i=1:nnodes_foot;
% %    for ip=1:nfoot;
% %     Dfft_gruppo((ip-1)*(nnodes_foot*6)+(i-1)*6+1,1)=Uffx(i,ip);
% %     Dfft_gruppo((ip-1)*(nnodes_foot*6)+(i-1)*6+2,1)=Uffy(i,ip);
% %     Dfft_gruppo((ip-1)*(nnodes_foot*6)+(i-1)*6+3,1)=Uffz(i,ip);
% %    end
% % end
% % 



%% 
% =========================================================================
% STIFFNESS MATRIX OF FOOTINGS
% =========================================================================

KKfoot=zeros(nfoot*(length(dx)*6+6),nfoot*(length(dx)*6+6));
Kpalo=zeros(length(dx)*6+6,length(dx)*6+6,nfoot);

for i=1:nfoot
    
    for ii=1:length(dx);  
    KBern3Delt =KBern3D_foot_TIM(Efoot,dfoot,bfoot,dx(ii),EGratio,coeff_TIM,ni_foot);%KBern3D(Ep,d,dx(ii),alfa(i),beta(i));    KBern3Delt =KBern3D_foot_TIM(Efoot,dfoot,bfoot,dx(ii),coeff_TIM);%KBern3D(Ep,d,dx(ii),alfa(i),beta(i));
%     KBern3Delt =KBern3D_foot(Efoot,dfoot,bfoot,dx(ii));%KBern3D(Ep,d,dx(ii),alfa(i),beta(i));    KBern3Delt =KBern3D_foot_TIM(Efoot,dfoot,bfoot,dx(ii),coeff_TIM);%KBern3D(Ep,d,dx(ii),alfa(i),beta(i));
        
    
    Kpalo((ii-1)*6+1:(ii-1)*6+12,(ii-1)*6+1:(ii-1)*6+12,i)=Kpalo((ii-1)*6+1:(ii-1)*6+12,(ii-1)*6+1:(ii-1)*6+12,i)+KBern3Delt;
    end
    
    KKfoot((i-1)*((length(dx)*6+6))+1:(i)*((length(dx)*6+6)),(i-1)*((length(dx)*6+6))+1:(i)*((length(dx)*6+6)))=...
    KKfoot((i-1)*((length(dx)*6+6))+1:(i)*((length(dx)*6+6)),(i-1)*((length(dx)*6+6))+1:(i)*((length(dx)*6+6)))+...
    Kpalo(:,:,i);

end


%% CONDENSED STIFFNESS MATRIX of the STRUCTURE
% DESCRIPTIVE TEXT

if length(X_foot_centr)>1
    
    % %%CONDENSED MATRIX OF BUILODING ON FIXED CONSTRAINTS
    [ RR_fixhead_foot,Reactions_qbeam_fixhead_foot ] = f_frame2D_cond_tie_pile_head_foot_v2( num_bay_foot,num_storey,spanx_frame,spanz_frame,Eframe,bc,dc,bb,db,bf,df,foundn,qz_beam );
    %  [ RR_fixhead ] = f_frame2D_cond_tie_pile_head_v1( num_bay,num_storey,spanx_frame,spanz_frame,Eframe,bc,dc,bb,db,bf,df,foundn );
    
    
    % % %%CONDENSED MATRIX OF BUILODING ON HINGED CONSTRAINTS
    [ RR_hinhead_foot,Reactions_qbeam_hinhead_foot ] = f_frame2D_cond_hinge_pile_head_foot_v2( num_bay_foot,num_storey,spanx_frame,spanz_frame,Eframe,bc,dc,bb,db,bf,df,foundn,qz_beam );
    % [ RR_hinhead ] = f_frame2D_cond_hinge_pile_head_v1( num_bay,num_storey,spanx_frame,spanz_frame,Eframe,bc,dc,bb,db,bf,df,foundn );
    
    % % %%% EQUIVALENT BEAM
    % % % DESCRIPTIVE TEXT
    % %
    % [ EI_mat,EA_mat ] = f_beam_eq_coeff( ht,X,X(1),X(length(X)),Eframe,bf,df,bc,dc,bb,db,spanz_frame,spanx_frame,foundn,num_storey );
    % %
    % % %tie_pile_head
    % [ RR_fixhead_eq_beam] = f_beam2D_cond_tie_pile_head_varEI_v1( X,zeros(1,npali),EI_mat,EA_mat );
    % %
    
    % The stiffness matrix computed by f_frame2D function has for column the
    % reaction forces at pile heads due to unit constrained displacement. This
    % condensed stiffness matrix need to be reorganized
    KKframe_fixed_foot=zeros(6*nfoot*nnodes_foot);
    KKframe_hinged_foot=zeros(6*nfoot*nnodes_foot);
    
    for ip=1:nfoot
        
        
        KKframe_fixed_foot(1+(nelementi_foot/2*6):(nelementi_foot+1)*6:nfoot*(nelementi_foot+1)*6,(ip-1)*6*(nelementi_foot+1)+1+(nelementi_foot/2*6))=RR_fixhead_foot(1:3:nfoot*3,(ip-1)*3+1);
        KKframe_fixed_foot(3+(nelementi_foot/2*6):(nelementi_foot+1)*6:nfoot*(nelementi_foot+1)*6,(ip-1)*6*(nelementi_foot+1)+1+(nelementi_foot/2*6))=RR_fixhead_foot(2:3:nfoot*3,(ip-1)*3+1);
        KKframe_fixed_foot(5+(nelementi_foot/2*6):(nelementi_foot+1)*6:nfoot*(nelementi_foot+1)*6,(ip-1)*6*(nelementi_foot+1)+1+(nelementi_foot/2*6))=RR_fixhead_foot(3:3:nfoot*3,(ip-1)*3+1);
        
        KKframe_fixed_foot(1+(nelementi_foot/2*6):(nelementi_foot+1)*6:nfoot*(nelementi_foot+1)*6,(ip-1)*6*(nelementi_foot+1)+3+(nelementi_foot/2*6))=RR_fixhead_foot(1:3:nfoot*3,(ip-1)*3+2);
        KKframe_fixed_foot(3+(nelementi_foot/2*6):(nelementi_foot+1)*6:nfoot*(nelementi_foot+1)*6,(ip-1)*6*(nelementi_foot+1)+3+(nelementi_foot/2*6))=RR_fixhead_foot(2:3:nfoot*3,(ip-1)*3+2);
        KKframe_fixed_foot(5+(nelementi_foot/2*6):(nelementi_foot+1)*6:nfoot*(nelementi_foot+1)*6,(ip-1)*6*(nelementi_foot+1)+3+(nelementi_foot/2*6))=RR_fixhead_foot(3:3:nfoot*3,(ip-1)*3+2);
        
        KKframe_fixed_foot(1+(nelementi_foot/2*6):(nelementi_foot+1)*6:nfoot*(nelementi_foot+1)*6,(ip-1)*6*(nelementi_foot+1)+5+(nelementi_foot/2*6))=RR_fixhead_foot(1:3:nfoot*3,(ip-1)*3+3);
        KKframe_fixed_foot(3+(nelementi_foot/2*6):(nelementi_foot+1)*6:nfoot*(nelementi_foot+1)*6,(ip-1)*6*(nelementi_foot+1)+5+(nelementi_foot/2*6))=RR_fixhead_foot(2:3:nfoot*3,(ip-1)*3+3);
        KKframe_fixed_foot(5+(nelementi_foot/2*6):(nelementi_foot+1)*6:nfoot*(nelementi_foot+1)*6,(ip-1)*6*(nelementi_foot+1)+5+(nelementi_foot/2*6))=RR_fixhead_foot(3:3:nfoot*3,(ip-1)*3+3);
    end
    
    for ip=1:nfoot
        
        KKframe_hinged_foot(1+(nelementi_foot/2*6):(nelementi_foot+1)*6:nfoot*(nelementi_foot+1)*6,(ip-1)*6*(nelementi_foot+1)+1+(nelementi_foot/2*6))=RR_hinhead_foot(1:2:nfoot*2,(ip-1)*2+1);
        KKframe_hinged_foot(3+(nelementi_foot/2*6):(nelementi_foot+1)*6:nfoot*(nelementi_foot+1)*6,(ip-1)*6*(nelementi_foot+1)+1+(nelementi_foot/2*6))=RR_hinhead_foot(2:2:nfoot*2,(ip-1)*2+1);
        
        KKframe_hinged_foot(1+(nelementi_foot/2*6):(nelementi_foot+1)*6:nfoot*(nelementi_foot+1)*6,(ip-1)*6*(nelementi_foot+1)+3+(nelementi_foot/2*6))=RR_hinhead_foot(1:2:nfoot*2,(ip-1)*2+2);
        KKframe_hinged_foot(3+(nelementi_foot/2*6):(nelementi_foot+1)*6:nfoot*(nelementi_foot+1)*6,(ip-1)*6*(nelementi_foot+1)+3+(nelementi_foot/2*6))=RR_hinhead_foot(2:2:nfoot*2,(ip-1)*2+2);
        
    end

    
   
elseif length(X_foot_centr)==1 && max(spanx_frame)>0&& max(spanz_frame)>0
      
    
   
    
    
     % %%CONDENSED MATRIX OF BUILODING ON FIXED CONSTRAINTS
	[ RR_fixhead_foot,Reactions_qbeam_fixhead_foot,FL_qbeam_fixhead_foot ,KG_el_qbeam_fixhead_foot,Element_nodes,xx,yy,EI_mat,EA_mat,U_qbeam_fixhead_foot] = f_frame2D_cond_tie_pile_head_foot_v3( num_bay_raft,num_storey,spanx_frame,spanz_frame,Eframe,bc,dc,bb,db,bf,df,foundn,qz_beam,pz_col );
% 	[ RR_fixhead_foot,Reactions_qbeam_fixhead_foot ] = f_frame2D_cond_tie_pile_head_foot_v2( num_bay_raft,num_storey,spanx_frame,spanz_frame,Eframe,bc,dc,bb,db,bf,df,foundn,qz_beam );
    %  [ RR_fixhead ] = f_frame2D_cond_tie_pile_head_v1( num_bay,num_storey,spanx_frame,spanz_frame,Eframe,bc,dc,bb,db,bf,df,foundn );
    
    
    % % %%CONDENSED MATRIX OF BUILODING ON HINGED CONSTRAINTS
    [ RR_hinhead_foot,Reactions_qbeam_hinhead_foot,FL_qbeam_hinhead_foot,KG_el_qbeam_hinhead_foot ] = f_frame2D_cond_hinge_pile_head_foot_v3( num_bay_raft,num_storey,spanx_frame,spanz_frame,Eframe,bc,dc,bb,db,bf,df,foundn,qz_beam,pz_col );
%     [ RR_hinhead_foot,Reactions_qbeam_hinhead_foot ] = f_frame2D_cond_hinge_pile_head_foot_v2( num_bay_raft,num_storey,spanx_frame,spanz_frame,Eframe,bc,dc,bb,db,bf,df,foundn,qz_beam );
    % [ RR_hinhead ] = f_frame2D_cond_hinge_pile_head_v1( num_bay,num_storey,spanx_frame,spanz_frame,Eframe,bc,dc,bb,db,bf,df,foundn );
    
    % % %%% EQUIVALENT BEAM
    % % % DESCRIPTIVE TEXT
    % %
    % [ EI_mat,EA_mat ] = f_beam_eq_coeff( ht,X,X(1),X(length(X)),Eframe,bf,df,bc,dc,bb,db,spanz_frame,spanx_frame,foundn,num_storey );
    % %
    % % %tie_pile_head
    % [ RR_fixhead_eq_beam] = f_beam2D_cond_tie_pile_head_varEI_v1( X,zeros(1,npali),EI_mat,EA_mat );
    % %
    
    % The stiffness matrix computed by f_frame2D function has for column the
    % reaction forces at pile heads due to unit constrained displacement. This
    % condensed stiffness matrix need to be reorganized
    KKframe_fixed_foot=zeros(6*nfoot*nnodes_foot);
    KKframe_hinged_foot=zeros(6*nfoot*nnodes_foot);
    
    for ip=1:(num_bay_raft+1)
        


        KKframe_fixed_foot(1:(nelementi_bayraft)*6:end,(ip-1)*6*(nelementi_bayraft)+1)=RR_fixhead_foot(1:3:(num_bay_raft+1)*3,(ip-1)*3+1);
        KKframe_fixed_foot(3:(nelementi_bayraft)*6:end,(ip-1)*6*(nelementi_bayraft)+1)=RR_fixhead_foot(2:3:(num_bay_raft+1)*3,(ip-1)*3+1);
        KKframe_fixed_foot(5:(nelementi_bayraft)*6:end,(ip-1)*6*(nelementi_bayraft)+1)=RR_fixhead_foot(3:3:(num_bay_raft+1)*3,(ip-1)*3+1);
        
        KKframe_fixed_foot(1:(nelementi_bayraft)*6:end,(ip-1)*6*(nelementi_bayraft)+3)=RR_fixhead_foot(1:3:(num_bay_raft+1)*3,(ip-1)*3+2);
        KKframe_fixed_foot(3:(nelementi_bayraft)*6:end,(ip-1)*6*(nelementi_bayraft)+3)=RR_fixhead_foot(2:3:(num_bay_raft+1)*3,(ip-1)*3+2);
        KKframe_fixed_foot(5:(nelementi_bayraft)*6:end,(ip-1)*6*(nelementi_bayraft)+3)=RR_fixhead_foot(3:3:(num_bay_raft+1)*3,(ip-1)*3+2);
        
        KKframe_fixed_foot(1:(nelementi_bayraft)*6:end,(ip-1)*6*(nelementi_bayraft)+5)=RR_fixhead_foot(1:3:(num_bay_raft+1)*3,(ip-1)*3+3);
        KKframe_fixed_foot(3:(nelementi_bayraft)*6:end,(ip-1)*6*(nelementi_bayraft)+5)=RR_fixhead_foot(2:3:(num_bay_raft+1)*3,(ip-1)*3+3);
        KKframe_fixed_foot(5:(nelementi_bayraft)*6:end,(ip-1)*6*(nelementi_bayraft)+5)=RR_fixhead_foot(3:3:(num_bay_raft+1)*3,(ip-1)*3+3);
    end
    
    for ip=1:(num_bay_raft+1)
        
        KKframe_hinged_foot(1:(nelementi_bayraft)*6:end,(ip-1)*6*(nelementi_bayraft)+1)=RR_hinhead_foot(1:2:(num_bay_raft+1)*2,(ip-1)*2+1);
        KKframe_hinged_foot(3:(nelementi_bayraft)*6:end,(ip-1)*6*(nelementi_bayraft)+1)=RR_hinhead_foot(2:2:(num_bay_raft+1)*2,(ip-1)*2+1);
        
        KKframe_hinged_foot(1:(nelementi_bayraft)*6:end,(ip-1)*6*(nelementi_bayraft)+3)=RR_hinhead_foot(1:2:(num_bay_raft+1)*2,(ip-1)*2+2);
        KKframe_hinged_foot(3:(nelementi_bayraft)*6:end,(ip-1)*6*(nelementi_bayraft)+3)=RR_hinhead_foot(2:2:(num_bay_raft+1)*2,(ip-1)*2+2);
        
    end 
     
    
    
elseif length(X_foot_centr)==1 && spanx_frame==0&& spanz_frame==0
    KKframe_fixed_foot=zeros(6*nfoot*nnodes_foot);
    KKframe_hinged_foot=zeros(6*nfoot*nnodes_foot);
    Reactions_qbeam_fixhead_foot=zeros(3*nfoot,1);
    Reactions_qbeam_hinhead_foot=zeros(3*nfoot,1);
    
else
    
msg = 'Error occurred in the definition of the frame.';
error(msg)
    
end




%%
% =========================================================================
% CALCOLO DELLA MATRICE DI VINCOLO A
% =========================================================================

%d=AdE
%d=[s1,...,sp,...,sn]t   dove sp spostamenti nodali del pali p-esimo
%dE[sE1,...,sEp,...,sEn]t   dove sEp spostamenti nodali del pali p-esimo a
%meno degli spostamenti a quota z=0

% A_cap_i=zeros((6*length(dz)+6)*npali,6+(6*length(dz))*npali);
% A_cap_h=zeros((6*length(dz)+6)*npali,6+(6*length(dz)+3)*npali);
% 
% for i=1:npali;
%     A_cap_i((i-1)*(6*length(dz)+6)+1:(i-1)*(6*length(dz)+6)+6,1:6)=Ai(Xm,Ym,Zm,X(i),Y(i),z(1));   
%     A_cap_i(7+(i-1)*(6*length(dz)+6):(i)*(6*length(dz)+6),7+(i-1)*(6*length(dz)):6+(i)*(6*length(dz)))=eye(6*length(dz),6*length(dz));
% end
% 
% for i=1:npali;
%     A_cap_h((i-1)*(6*length(dz)+6)+1:(i-1)*(6*length(dz)+6)+3,1:6)=Ah(Xm,Ym,Zm,X(i),Y(i),z(1));   
%     A_cap_h(4+(i-1)*(6*length(dz)+6):(i)*(6*length(dz)+6),7+(i-1)*(6*length(dz)+3):6+(i)*(6*length(dz)+3))=eye(6*length(dz)+3,6*length(dz)+3);
% end


%% 
% =========================================================================
% FLEXIBILITY MATRIX COMPUTED WITH VAZIRI ET AL. (1982)
% =========================================================================

FLEX=zeros(nfoot*(nelementi_foot+1)*3);


%FLEXIBILITY MATRIX
for jg=1:nfoot*(nelementi_foot+1)
    
    
    Zj_nod(1:1:nfoot*(nelementi_foot+1),1)=Zi_nod(jg);
    Xj_nod(1:1:nfoot*(nelementi_foot+1),1)=Xi_nod(jg);
    Yj_nod(1:1:nfoot*(nelementi_foot+1),1)=Yi_nod(jg);

    dZj_nod(:,jg)=Zi_nod(:)-Zj_nod(:);
    dXj_nod(:,jg)=Xi_nod(:)-Xj_nod(:);
    dYj_nod(:,jg)=Yi_nod(:)-Yj_nod(:); 
    
    [ uxij_1,uyij_1,uzij_1 ] = int_factor_mindlin_j_1_vect_Vasiri( Xi_nod,Yi_nod,Zi_nod,Xj_nod,Yj_nod,Zj_nod,Gs_nod,ni_nod,h_el_foot,bfoot);
    [ uxij_2,uyij_2,uzij_2 ] = int_factor_mindlin_j_2_vect_Vasiri( Xi_nod,Yi_nod,Zi_nod,Xj_nod,Yj_nod,Zj_nod,Gs_nod,ni_nod,h_el_foot,bfoot);
    [ uxij_3,uyij_3,uzij_3 ] = int_factor_mindlin_j_3_vect_Vasiri( Xi_nod,Yi_nod,Zi_nod,Xj_nod,Yj_nod,Zj_nod,Gs_nod,ni_nod,h_el_foot,bfoot);

    FLEX(1:3:nfoot*(nelementi_foot+1)*3,(jg-1)*3+1)=uxij_1;
    FLEX(2:3:nfoot*(nelementi_foot+1)*3,(jg-1)*3+1)=uyij_1;
    FLEX(3:3:nfoot*(nelementi_foot+1)*3,(jg-1)*3+1)=uzij_1;
    
    FLEX(1:3:nfoot*(nelementi_foot+1)*3,(jg-1)*3+2)=uxij_2;
    FLEX(2:3:nfoot*(nelementi_foot+1)*3,(jg-1)*3+2)=uyij_2;
    FLEX(3:3:nfoot*(nelementi_foot+1)*3,(jg-1)*3+2)=uzij_2;

    FLEX(1:3:nfoot*(nelementi_foot+1)*3,(jg-1)*3+3)=uxij_3;
    FLEX(2:3:nfoot*(nelementi_foot+1)*3,(jg-1)*3+3)=uyij_3;
    FLEX(3:3:nfoot*(nelementi_foot+1)*3,(jg-1)*3+3)=uzij_3;

end

%Coefficients to consider that the influence aread at the edges of the
%footings is half
r1 = -0.3136*exp(-0.531*(h_el_foot/bfoot))+0.2174*(h_el_foot/bfoot)^0.1755+1.327;
r2 = -0.3136*exp(-0.531*(bfoot/h_el_foot))+0.2174*(bfoot/h_el_foot)^0.1755+1.327;
r3 = -0.3136*exp(-0.531*(h_el_foot/bfoot))+0.2174*(h_el_foot/bfoot)^0.1755+1.327;

for ip=1:nfoot
    FLEX(1+(ip-1)*nnodes_foot*3,1+(ip-1)*nnodes_foot*3)=FLEX(1+(ip-1)*nnodes_foot*3,1+(ip-1)*nnodes_foot*3)*r1;
    FLEX(2+(ip-1)*nnodes_foot*3,2+(ip-1)*nnodes_foot*3)=FLEX(2+(ip-1)*nnodes_foot*3,2+(ip-1)*nnodes_foot*3)*r2;
    FLEX(3+(ip-1)*nnodes_foot*3,3+(ip-1)*nnodes_foot*3)=FLEX(3+(ip-1)*nnodes_foot*3,3+(ip-1)*nnodes_foot*3)*r3;
    FLEX(1+(nnodes_foot-1)*3+(ip-1)*nnodes_foot*3,1+(nnodes_foot-1)*3+(ip-1)*nnodes_foot*3)=FLEX(1+(nnodes_foot-1)*3+(ip-1)*nnodes_foot*3,1+(nnodes_foot-1)*3+(ip-1)*nnodes_foot*3)*r1;
    FLEX(2+(nnodes_foot-1)*3+(ip-1)*nnodes_foot*3,2+(nnodes_foot-1)*3+(ip-1)*nnodes_foot*3)=FLEX(2+(nnodes_foot-1)*3+(ip-1)*nnodes_foot*3,2+(nnodes_foot-1)*3+(ip-1)*nnodes_foot*3)*r2;
    FLEX(3+(nnodes_foot-1)*3+(ip-1)*nnodes_foot*3,3+(nnodes_foot-1)*3+(ip-1)*nnodes_foot*3)=FLEX(3+(nnodes_foot-1)*3+(ip-1)*nnodes_foot*3,3+(nnodes_foot-1)*3+(ip-1)*nnodes_foot*3)*r3;
end

% for ip=1:nfoot
%     FLEX(1+(ip-1)*nnodes_foot*3,1+(ip-1)*nnodes_foot*3)=FLEX(1+(ip-1)*nnodes_foot*3,1+(ip-1)*nnodes_foot*3)*2;
%     FLEX(2+(ip-1)*nnodes_foot*3,2+(ip-1)*nnodes_foot*3)=FLEX(2+(ip-1)*nnodes_foot*3,2+(ip-1)*nnodes_foot*3)*2;
%     FLEX(3+(ip-1)*nnodes_foot*3,3+(ip-1)*nnodes_foot*3)=FLEX(3+(ip-1)*nnodes_foot*3,3+(ip-1)*nnodes_foot*3)*2;
%     FLEX(1+(nnodes_foot-1)*3+(ip-1)*nnodes_foot*3,1+(nnodes_foot-1)*3+(ip-1)*nnodes_foot*3)=FLEX(1+(nnodes_foot-1)*3+(ip-1)*nnodes_foot*3,1+(nnodes_foot-1)*3+(ip-1)*nnodes_foot*3)*2;
%     FLEX(2+(nnodes_foot-1)*3+(ip-1)*nnodes_foot*3,2+(nnodes_foot-1)*3+(ip-1)*nnodes_foot*3)=FLEX(2+(nnodes_foot-1)*3+(ip-1)*nnodes_foot*3,2+(nnodes_foot-1)*3+(ip-1)*nnodes_foot*3)*2;
%     FLEX(3+(nnodes_foot-1)*3+(ip-1)*nnodes_foot*3,3+(nnodes_foot-1)*3+(ip-1)*nnodes_foot*3)=FLEX(3+(nnodes_foot-1)*3+(ip-1)*nnodes_foot*3,3+(nnodes_foot-1)*3+(ip-1)*nnodes_foot*3)*2;
% end

el_area=bfoot*h_el_foot*ones(nfoot*nnodes_foot,1);
for ip=1:nfoot
    el_area(1+(ip-1)*nnodes_foot)=0.5*el_area(1+(ip-1)*nnodes_foot);
    el_area(1+(nnodes_foot-1)+(ip-1)*nnodes_foot)=0.5*el_area(1+(nnodes_foot-1)+(ip-1)*nnodes_foot);
end

 %% 
% % =========================================================================
% % FLEXIBILITY MATRIX
% %	GAZETAS (1991) FOR DIAGNOAL LOCAL STIFFNESS TERMS (assuming uniform
% %	displacement beneath finite elements)
% %   MINDLIND(1931) FOR OFF DIAGONAL INTERACTION TERMS
% % =========================================================================
% 
% FLEX_GAZ=zeros(nfoot*(nelementi_foot+1)*3);
% 
% %Winkler model of the soil
[ kh_gazetas,kv_gazetas ] = soilsprings_static_foot_gazetas(Es(1),Efoot,nis,bfoot,dfoot,h_el_foot); 
% 
% %FLEXIBILITY MATRIX
% for jg=1:nfoot*(nelementi_foot+1)
%     
%     Zj_nod(1:1:nfoot*(nelementi_foot+1),1)=Zi_nod(jg);
%     Xj_nod(1:1:nfoot*(nelementi_foot+1),1)=Xi_nod(jg);
%     Yj_nod(1:1:nfoot*(nelementi_foot+1),1)=Yi_nod(jg);
% 
%     dZj_nod(:,jg)=Zi_nod(:)-Zj_nod(:);
%     dXj_nod(:,jg)=Xi_nod(:)-Xj_nod(:);
%     dYj_nod(:,jg)=Yi_nod(:)-Yj_nod(:); 
%    
%     [ uxij_1,uyij_1,uzij_1 ] = int_factor_mindlin_j_1_vect_gaz( Xi_nod,Yi_nod,Zi_nod,Xj_nod,Yj_nod,Zj_nod,Gs_nod,ni_nod,kh_gazetas,kv_gazetas,h_el_foot,nfoot,nnodes_foot);
%     [ uxij_2,uyij_2,uzij_2 ] = int_factor_mindlin_j_2_vect_gaz( Xi_nod,Yi_nod,Zi_nod,Xj_nod,Yj_nod,Zj_nod,Gs_nod,ni_nod,kh_gazetas,kv_gazetas,h_el_foot,nfoot,nnodes_foot);
%     [ uxij_3,uyij_3,uzij_3 ] = int_factor_mindlin_j_3_vect_gaz( Xi_nod,Yi_nod,Zi_nod,Xj_nod,Yj_nod,Zj_nod,Gs_nod,ni_nod,kh_gazetas,kv_gazetas,h_el_foot,nfoot,nnodes_foot);
% 
%     FLEX_GAZ(1:3:nfoot*(nelementi_foot+1)*3,(jg-1)*3+1)=uxij_1;
%     FLEX_GAZ(2:3:nfoot*(nelementi_foot+1)*3,(jg-1)*3+1)=uyij_1;
%     FLEX_GAZ(3:3:nfoot*(nelementi_foot+1)*3,(jg-1)*3+1)=uzij_1;
%     
%     FLEX_GAZ(1:3:nfoot*(nelementi_foot+1)*3,(jg-1)*3+2)=uxij_2;
%     FLEX_GAZ(2:3:nfoot*(nelementi_foot+1)*3,(jg-1)*3+2)=uyij_2;
%     FLEX_GAZ(3:3:nfoot*(nelementi_foot+1)*3,(jg-1)*3+2)=uzij_2;
% 
%     FLEX_GAZ(1:3:nfoot*(nelementi_foot+1)*3,(jg-1)*3+3)=uxij_3;
%     FLEX_GAZ(2:3:nfoot*(nelementi_foot+1)*3,(jg-1)*3+3)=uyij_3;
%     FLEX_GAZ(3:3:nfoot*(nelementi_foot+1)*3,(jg-1)*3+3)=uzij_3;
% 
% end
% 
% for ip=1:nfoot
%     FLEX_GAZ(1+(ip-1)*nnodes_foot*3,1+(ip-1)*nnodes_foot*3)=FLEX_GAZ(1+(ip-1)*nnodes_foot*3,1+(ip-1)*nnodes_foot*3)*2;
%     FLEX_GAZ(2+(ip-1)*nnodes_foot*3,2+(ip-1)*nnodes_foot*3)=FLEX_GAZ(2+(ip-1)*nnodes_foot*3,2+(ip-1)*nnodes_foot*3)*2;
%     FLEX_GAZ(3+(ip-1)*nnodes_foot*3,3+(ip-1)*nnodes_foot*3)=FLEX_GAZ(3+(ip-1)*nnodes_foot*3,3+(ip-1)*nnodes_foot*3)*2;
%     FLEX_GAZ(1+(nnodes_foot-1)*3+(ip-1)*nnodes_foot*3,1+(nnodes_foot-1)*3+(ip-1)*nnodes_foot*3)=FLEX_GAZ(1+(nnodes_foot-1)*3+(ip-1)*nnodes_foot*3,1+(nnodes_foot-1)*3+(ip-1)*nnodes_foot*3)*2;
%     FLEX_GAZ(2+(nnodes_foot-1)*3+(ip-1)*nnodes_foot*3,2+(nnodes_foot-1)*3+(ip-1)*nnodes_foot*3)=FLEX_GAZ(2+(nnodes_foot-1)*3+(ip-1)*nnodes_foot*3,2+(nnodes_foot-1)*3+(ip-1)*nnodes_foot*3)*2;
%     FLEX_GAZ(3+(nnodes_foot-1)*3+(ip-1)*nnodes_foot*3,3+(nnodes_foot-1)*3+(ip-1)*nnodes_foot*3)=FLEX_GAZ(3+(nnodes_foot-1)*3+(ip-1)*nnodes_foot*3,3+(nnodes_foot-1)*3+(ip-1)*nnodes_foot*3)*2;
% end

%% 
% =========================================================================
% FLEXIBILITY MATRIX COMPUTED WITH VESIC (1991) AND NO INTERACTION OFF DIAGONAL TERMS
% =========================================================================

FLEX_Vesic=zeros(nfoot*(nelementi_foot+1)*3);

%Winkler model of the soil
[ ksh_vesic,ksv_vesic ] = soilsprings_static_foot_vesic(Es(1),Efoot,nis,bfoot,dfoot); 

%FLEXIBILITY MATRIX
for jg=1:nfoot*(nelementi_foot+1)

    Zj_nod(1:1:nfoot*(nelementi_foot+1),1)=Zi_nod(jg);
    Xj_nod(1:1:nfoot*(nelementi_foot+1),1)=Xi_nod(jg);
    Yj_nod(1:1:nfoot*(nelementi_foot+1),1)=Yi_nod(jg);

    dZj_nod(:,jg)=Zi_nod(:)-Zj_nod(:);
    dXj_nod(:,jg)=Xi_nod(:)-Xj_nod(:);
    dYj_nod(:,jg)=Yi_nod(:)-Yj_nod(:); 

    [ uxij_1,uyij_1,uzij_1 ] = int_factor_mindlin_j_1_vect_jap( Xi_nod,Yi_nod,Zi_nod,Xj_nod,Yj_nod,Zj_nod,Gs_nod,ni_nod,ksh_vesic,ksv_vesic,h_el_foot,nfoot,nnodes_foot);
    [ uxij_2,uyij_2,uzij_2 ] = int_factor_mindlin_j_2_vect_jap( Xi_nod,Yi_nod,Zi_nod,Xj_nod,Yj_nod,Zj_nod,Gs_nod,ni_nod,ksh_vesic,ksv_vesic,h_el_foot,nfoot,nnodes_foot);
    [ uxij_3,uyij_3,uzij_3 ] = int_factor_mindlin_j_3_vect_jap( Xi_nod,Yi_nod,Zi_nod,Xj_nod,Yj_nod,Zj_nod,Gs_nod,ni_nod,ksh_vesic,ksv_vesic,h_el_foot,nfoot,nnodes_foot);
    
    FLEX_Vesic(1:3:nfoot*(nelementi_foot+1)*3,(jg-1)*3+1)=uxij_1;
    FLEX_Vesic(2:3:nfoot*(nelementi_foot+1)*3,(jg-1)*3+1)=uyij_1;
    FLEX_Vesic(3:3:nfoot*(nelementi_foot+1)*3,(jg-1)*3+1)=uzij_1;
    
    FLEX_Vesic(1:3:nfoot*(nelementi_foot+1)*3,(jg-1)*3+2)=uxij_2;
    FLEX_Vesic(2:3:nfoot*(nelementi_foot+1)*3,(jg-1)*3+2)=uyij_2;
    FLEX_Vesic(3:3:nfoot*(nelementi_foot+1)*3,(jg-1)*3+2)=uzij_2;

    FLEX_Vesic(1:3:nfoot*(nelementi_foot+1)*3,(jg-1)*3+3)=uxij_3;
    FLEX_Vesic(2:3:nfoot*(nelementi_foot+1)*3,(jg-1)*3+3)=uyij_3;
    FLEX_Vesic(3:3:nfoot*(nelementi_foot+1)*3,(jg-1)*3+3)=uzij_3;

end

% for ip=1:nfoot
%     FLEX_Vesic(1+(ip-1)*nnodes_foot*3,1+(ip-1)*nnodes_foot*3)=FLEX_Vesic(1+(ip-1)*nnodes_foot*3,1+(ip-1)*nnodes_foot*3)*2;
%     FLEX_Vesic(2+(ip-1)*nnodes_foot*3,2+(ip-1)*nnodes_foot*3)=FLEX_Vesic(2+(ip-1)*nnodes_foot*3,2+(ip-1)*nnodes_foot*3)*2;
%     FLEX_Vesic(3+(ip-1)*nnodes_foot*3,3+(ip-1)*nnodes_foot*3)=FLEX_Vesic(3+(ip-1)*nnodes_foot*3,3+(ip-1)*nnodes_foot*3)*2;
%     FLEX_Vesic(1+(nnodes_foot-1)*3+(ip-1)*nnodes_foot*3,1+(nnodes_foot-1)*3+(ip-1)*nnodes_foot*3)=FLEX_Vesic(1+(nnodes_foot-1)*3+(ip-1)*nnodes_foot*3,1+(nnodes_foot-1)*3+(ip-1)*nnodes_foot*3)*2;
%     FLEX_Vesic(2+(nnodes_foot-1)*3+(ip-1)*nnodes_foot*3,2+(nnodes_foot-1)*3+(ip-1)*nnodes_foot*3)=FLEX_Vesic(2+(nnodes_foot-1)*3+(ip-1)*nnodes_foot*3,2+(nnodes_foot-1)*3+(ip-1)*nnodes_foot*3)*2;
%     FLEX_Vesic(3+(nnodes_foot-1)*3+(ip-1)*nnodes_foot*3,3+(nnodes_foot-1)*3+(ip-1)*nnodes_foot*3)=FLEX_Vesic(3+(nnodes_foot-1)*3+(ip-1)*nnodes_foot*3,3+(nnodes_foot-1)*3+(ip-1)*nnodes_foot*3)*2;
% end

% FLEXIBILITY MATRIX WITHOUT INTERACTION BETWEEN NODES
FLEX_Vesic=diag(diag(FLEX_Vesic));

% FLEXIBILITY MATRIX WITHOUT INTERACTION BETWEEN NODES IN THE SAME PILE
% Interection effects between the nodes within the same pile is neglected:
% components outside diagonal are set to zero if nodes i and j with same
% pile

% for in=1:(nfoot*(nelementi_foot+1)*3)
%     for jn=1:(nfoot*(nelementi_foot+1)*3)
%         i=ceil(in/(nelementi_foot+1)/3);  %ie
%         j=ceil(jn/(nelementi_foot+1)/3);  %je
%         if i==j &&  in~=jn
%         FLEX_Vesic(in,jn)=0;
%         end
%     end
% end


for ip=1:nfoot
    FLEX_Vesic(1+(ip-1)*nnodes_foot*3,1+(ip-1)*nnodes_foot*3)=FLEX_Vesic(1+(ip-1)*nnodes_foot*3,1+(ip-1)*nnodes_foot*3)*2;
    FLEX_Vesic(2+(ip-1)*nnodes_foot*3,2+(ip-1)*nnodes_foot*3)=FLEX_Vesic(2+(ip-1)*nnodes_foot*3,2+(ip-1)*nnodes_foot*3)*2;
    FLEX_Vesic(3+(ip-1)*nnodes_foot*3,3+(ip-1)*nnodes_foot*3)=FLEX_Vesic(3+(ip-1)*nnodes_foot*3,3+(ip-1)*nnodes_foot*3)*2;
    FLEX_Vesic(1+(nnodes_foot-1)*3+(ip-1)*nnodes_foot*3,1+(nnodes_foot-1)*3+(ip-1)*nnodes_foot*3)=FLEX_Vesic(1+(nnodes_foot-1)*3+(ip-1)*nnodes_foot*3,1+(nnodes_foot-1)*3+(ip-1)*nnodes_foot*3)*2;
    FLEX_Vesic(2+(nnodes_foot-1)*3+(ip-1)*nnodes_foot*3,2+(nnodes_foot-1)*3+(ip-1)*nnodes_foot*3)=FLEX_Vesic(2+(nnodes_foot-1)*3+(ip-1)*nnodes_foot*3,2+(nnodes_foot-1)*3+(ip-1)*nnodes_foot*3)*2;
    FLEX_Vesic(3+(nnodes_foot-1)*3+(ip-1)*nnodes_foot*3,3+(nnodes_foot-1)*3+(ip-1)*nnodes_foot*3)=FLEX_Vesic(3+(nnodes_foot-1)*3+(ip-1)*nnodes_foot*3,3+(nnodes_foot-1)*3+(ip-1)*nnodes_foot*3)*2;
end


%% Debug configuration
% %  FLEX=FLEX_Vesic;
 


%%
%STIFFNESS MATRIX OF THE SOIL
%Only translational dofs
KKsoil_rid=inv(FLEX);

if exist('switch_Ritterfoundation')==1
if switch_Ritterfoundation==1

Vn=1;
Hn=1;
ni=nis;
a=75*10.7/1000;%b_xloc(jg);
b=75*78.6/2/1000;%b_yloc(jg);
n=0;
m=0;
B=(2*(m-n)+1).*a./b;
C=(2*(m-n)-1).*a./b;
wmn=2.*Vn.*(1-ni.^2)/a/pi/Es.*(B.*asinh(1./B)+asinh(B));
uxMN=2.*Hn.*(1-ni.^2)./pi./Es./a.*(B.*asinh(1./B)+asinh(B))+...
    2.*Hn.*ni.*(1+ni)./pi./Es./a.*(asinh(B));
   
    
KKsoil_rid(3,3)=KKsoil_rid(3,3)+1/wmn;
KKsoil_rid(end,end)=KKsoil_rid(end,end)+1/wmn;

KKsoil_rid(1,1)=KKsoil_rid(1,1)+1/uxMN;
KKsoil_rid(end-2,end-2)=KKsoil_rid(end-2,end-2)+1/uxMN;

FLEX=inv(KKsoil_rid);

end
end



KKsoil=zeros(nfoot*(nelementi_foot+1)*6);
%Zeros are added within soil stiffness matrix to account for rotational
%dofs
for i=1:nfoot*(nelementi_foot+1)
    for j=1:nfoot*(nelementi_foot+1)
       KKsoil(1+(i-1)*6:3+(i-1)*6,1+(j-1)*6:3+(j-1)*6)=...
        KKsoil_rid(1+(i-1)*3:3+(i-1)*3,1+(j-1)*3:3+(j-1)*3);
	end
end

%ADDITIONAL SPRINGS IN THE DIRECITONS 2 AND 4 TO REMOVE UNCOSTRAINED DEGREE OF FREEDOM 
KKextra=zeros(6*nnodes_foot*nfoot);
KKextra(4:6:6*nnodes_foot*nfoot,4:6:6*nnodes_foot*nfoot)=kv_gazetas*10^3;
KKextra(2:6:6*nnodes_foot*nfoot,2:6:6*nnodes_foot*nfoot)=kv_gazetas*10^3;
KKextra=diag(diag(KKextra));

%GLOBAL STIFFNESS MATRIX    
    KKpg=KKfoot+KKsoil+KKextra;  
    
%GLOBAL STIFFNESS MATRIX FRAME WITH FIXED PILE HEADS  
    KKpg_fr_i=KKframe_fixed_foot+KKfoot+KKsoil+KKextra; 

%GLOBAL STIFFNESS MATRIX FRAME WITH HINGED PILE HEADS  
    KKpg_fr_h=KKframe_hinged_foot+KKfoot+KKsoil+KKextra; 
    
% % %EXTERNAL FORCES
% % Fpg=(KKsoil)*Dfft_gruppo; %vettore forze esterne


%% 
% =========================================================================
% INPUTS
% =========================================================================


% ------------------------------------------------------------------------- 
% PARAMETERS 
% ------------------------------------------------------------------------- 



    num_incr_TUNNEL=20;
    num_incr_P=num_incr_TUNNEL/2; 
    
    beta=1;
    delta_err   =1.0e-3; %set at least at 1.0e-3
    eps_err     =1.0e-3; %set at least at 1.0e-3
    flim2       =0;         %maximum tensile force of plastic elements (<=0)
    flim1       =Inf;%It needs to be set as a limit prssure rather than a force considering area corresponding to each node%2.5*10^5/(0.5*1.2);%       %maximum compressive force of plastic elements   (>=0) %     flim1       = c *Nc *sc +  0.5 gamma * B *Ng *sg;
    mu=tan(2*pi/360*phi_int);    %friction coefficient
   
% ------------------------------------------------------------------------- 
% GREENFIELD TUNNELLING MOVEMENTS
% -------------------------------------------------------------------------   

if switch_tun_adv==1

 
%flag=0;  % Flag=0 allows plotting only critical tunnel drive plots

   
% Tunnelling parameters
% % 	flagclay=1; %flagclay=1 changes x and y greenfield for clay conditions
yf=-1000; %Distance of origin to tunnel portal
y0=+icdf(makedist('Normal',0,Ky*ht),delta_ratio);% Offset distance

% Building parameters  
theta=theta_degree/180*pi; % 36 degrees for  P2-P3 facade, 26 degrees camos
    
num_incr_tunadv=length(yslim1:ysstep:yslim2);
ys_vect=yslim1:ysstep:yslim2;
if num_incr_tunadv<10 
h = msgbox('The minimum number of tunnel advancements is 10 for the elastoplastic solution', 'Error','error');
stop
end

%lines below are needed so that the following code recognises each adv as a tunnel vlt increment
num_incr_TUNNEL=num_incr_tunadv;
Vltp_vect=Vltp*ones(1,num_incr_tunadv+1);

ii=1;
for ys=yslim1:ysstep:yslim2 %Distance of origin to tunnel head
    % Building movements
    
    % Building coords local/global
    if nfoot==1
        lbuild=(max(x_global_foot(end))-min(x_global_foot(1)));
        xbl=x_global_foot-(x_global_foot(end)+x_global_foot(1))/2;
        x_global_foot=xbl+dorig+(lbuild/2); %Foundation coordinate in the xr and yr reference system (equal to the global rotated by theta)
        xnod_tot=xbl(:)+dorig+(lbuild/2);   %Foundation coordinate in the xr and yr reference system (equal to the global rotated by theta)
    else
        lbuild=(max(x_global_foot(end,end))-min(x_global_foot(1,1)));
        xbl=x_global_foot-(x_global_foot(end,end)+x_global_foot(1,1))/2;
        x_global_foot=xbl+dorig+(lbuild/2); %Foundation coordinate in the xr and yr reference system (equal to the global rotated by theta)
        xnod_tot=xbl(:)+dorig+(lbuild/2);   %Foundation coordinate in the xr and yr reference system (equal to the global rotated by theta)
    end
    xbg=forig+(dorig+(lbuild/2+xbl)).*cos(theta); % Global coord for grid
    ybg=(dorig+(lbuild/2+xbl)).*sin(theta);
    
    % Compute the 3D trough
    % For building
    [dxbg,dybg,dzbg]=comp3_AF(xbg,ybg,0,ys,yf,y0,Kx,Ky,ht,Vltp/100,Rt*2);
    dxbl=+cos(theta)*dxbg+sin(theta)*dybg;
    dybl=-sin(theta)*dxbg+cos(theta)*dybg;
    dzbl=dzbg;
    
    % Organise input displacement vector
    Uffx(:,:,ii)=dxbl;
    Uffz(:,:,ii)=dzbl;
    %     Uffy(:,:,ii)=dybl;%Note that Y movements could be added
    
    Uffx_vect(:,:,ii)=coeffux*Uffx(:,:,ii);
    Uffz_vect(:,:,ii)=coeffuz*Uffz(:,:,ii);
    Uffy_vect(:,:,ii)=zeros(nnodes_foot,nfoot);
    
    Uffz(isnan(Uffz)) = 0 ;
    Uffx(isnan(Uffx)) = 0 ;
    Uffx_vect(isnan(Uffx_vect)) = 0 ;
    Uffz_vect(isnan(Uffz_vect)) = 0 ;
    
    
    ii=ii+1;
end
I_gcurve_vect = 1000^2*Kx*ht*ones(1,num_incr_tunadv+1);



%Uff VECTOR AT NODES OF FINITE ELEMENTS
Dfft_gruppo_vect=zeros(nfoot*nnodes_foot*6,num_incr_tunadv);
for ii=1:length(yslim1:ysstep:yslim2)
    for i=1:nnodes_foot;
        for ip=1:nfoot;
            Dfft_gruppo_vect((ip-1)*(nnodes_foot*6)+(i-1)*6+1,ii)=Uffx_vect(i,ip,ii);
            Dfft_gruppo_vect((ip-1)*(nnodes_foot*6)+(i-1)*6+2,ii)=Uffy_vect(i,ip,ii);
            Dfft_gruppo_vect((ip-1)*(nnodes_foot*6)+(i-1)*6+3,ii)=Uffz_vect(i,ip,ii);
        end
    end
end

Dfft_gruppo_vect=[zeros(size(Dfft_gruppo_vect,1),1),Dfft_gruppo_vect];

vect_Ugf_z=Dfft_gruppo_vect(3:6:end,:);
vect_Ugf_x=Dfft_gruppo_vect(1:6:end,:);

            %*************TWIN TUNNEL
            if switch_tun_adv_twin==1
   
                
                %flag=0;  % Flag=0 allows plotting only critical tunnel drive plots
                
                
                % Tunnelling parameters
                % % 	flagclay=1; %flagclay=1 changes x and y greenfield for clay conditions
                yf_twin=-1000; %Distance of origin to tunnel portal
                y0_twin=+icdf(makedist('Normal',0,Ky_twin*ht_twin),delta_ratio_twin);% Offset distance
                
                % Building parameters
                theta_twin=theta_degree_twin/180*pi; % 36 degrees for  P2-P3 facade, 26 degrees camos
                
                num_incr_tunadv_twin=length(yslim1_twin:ysstep_twin:yslim2_twin);
                ys_vect_twin=yslim1_twin:ysstep_twin:yslim2_twin;
                if num_incr_tunadv_twin<10
                    h = msgbox('The minimum number of tunnel advancements is 10 for the elastoplastic solution', 'Error','error');
                    stop
                end
                
                %lines below are needed so that the following code recognises each adv as a tunnel vlt increment
                num_incr_TUNNEL_twin=num_incr_tunadv_twin;
                Vltp_vect_twin=Vltp_twin*ones(1,num_incr_tunadv_twin+1);
                
                ii=1;
                for ys_twin=yslim1_twin+ysstep_twin:ysstep_twin:yslim2_twin %Distance of origin to tunnel head
                    % Building movements
                    
                    % Building coords local/global
                    if nfoot==1
                        lbuild=(max(x_global_foot(end))-min(x_global_foot(1)));
                        xbl=x_global_foot-(x_global_foot(end)+x_global_foot(1))/2;
                        x_global_foot_twin=xbl+dorig_twin+(lbuild/2); %Foundation coordinate in the xr and yr reference system (equal to the global rotated by theta_twin)
                        xnod_tot_twin=xbl(:)+dorig_twin+(lbuild/2);   %Foundation coordinate in the xr and yr reference system (equal to the global rotated by theta_twin)
                    else
                        lbuild=(max(x_global_foot(end,end))-min(x_global_foot(1,1)));
                        xbl=x_global_foot-(x_global_foot(end,end)+x_global_foot(1,1))/2;
                        x_global_foot_twin=xbl+dorig_twin+(lbuild/2); %Foundation coordinate in the xr and yr reference system (equal to the global rotated by theta_twin)
                        xnod_tot_twin=xbl(:)+dorig_twin+(lbuild/2);   %Foundation coordinate in the xr and yr reference system (equal to the global rotated by theta_twin)
                    end
                    xbg_twin=forig_twin+(dorig_twin+(lbuild/2+xbl)).*cos(theta_twin); % Global coord for grid
                    ybg_twin=(dorig_twin+(lbuild/2+xbl)).*sin(theta_twin);
                    
                    % Compute the 3D trough
                    % For building
                    [dxbg_twin,dybg_twin,dzbg_twin]=comp3_AF(xbg_twin,ybg_twin,0,ys_twin,yf_twin,y0_twin,Kx_twin,Ky_twin,ht_twin,Vltp_twin/100,Rt_twin*2);
                    dxbl_twin=+cos(theta_twin)*dxbg_twin+sin(theta_twin)*dybg_twin;
                    dybl_twin=-sin(theta_twin)*dxbg_twin+cos(theta_twin)*dybg_twin;
                    dzbl_twin=dzbg_twin;
                    
                    % Organise input displacement vector
                    Uffx_twin(:,:,ii)=dxbl_twin;
                    Uffz_twin(:,:,ii)=dzbl_twin;
                    %     Uffy(:,:,ii)=dybl;%Note that Y movements could be added
                    
                    Uffx_vect_twin(:,:,ii)=coeffux_twin*Uffx_twin(:,:,ii);
                    Uffz_vect_twin(:,:,ii)=coeffuz_twin*Uffz_twin(:,:,ii);
                    Uffy_vect_twin(:,:,ii)=zeros(nnodes_foot,nfoot);
                    
                    Uffz_twin(isnan(Uffz_twin)) = 0 ;
                    Uffx_twin(isnan(Uffx_twin)) = 0 ;
                    Uffx_vect_twin(isnan(Uffx_vect_twin)) = 0 ;
                    Uffz_vect_twin(isnan(Uffz_vect_twin)) = 0 ;
                    
                    
                    ii=ii+1;
                end
                I_gcurve_vect_twin = 1000^2*Kx_twin*ht_twin*ones(1,num_incr_tunadv_twin+1);
                
                %Uff VECTOR AT NODES OF FINITE ELEMENTS
                Dfft_gruppo_vect_twin=zeros(nfoot*nnodes_foot*6,num_incr_tunadv_twin);
                for ii=1:length(yslim1_twin+ysstep_twin:ysstep_twin:yslim2_twin)
                    for i=1:nnodes_foot;
                        for ip=1:nfoot;
                            Dfft_gruppo_vect_twin((ip-1)*(nnodes_foot*6)+(i-1)*6+1,ii)=Uffx_vect_twin(i,ip,ii);
                            Dfft_gruppo_vect_twin((ip-1)*(nnodes_foot*6)+(i-1)*6+2,ii)=Uffy_vect_twin(i,ip,ii);
                            Dfft_gruppo_vect_twin((ip-1)*(nnodes_foot*6)+(i-1)*6+3,ii)=Uffz_vect_twin(i,ip,ii);
                        end
                    end
                end
                
                Dfft_gruppo_vect_twin=[zeros(size(Dfft_gruppo_vect_twin,1),1),Dfft_gruppo_vect_twin];
                
                vect_Ugf_z_twin=Dfft_gruppo_vect_twin(3:6:end,:);
                vect_Ugf_x_twin=Dfft_gruppo_vect_twin(1:6:end,:);
                
                
                
                
                
                vect_Ugf_x=[vect_Ugf_x,repmat(vect_Ugf_x(:,end),[1 size(vect_Ugf_x_twin,2)])+vect_Ugf_x_twin];
                vect_Ugf_z=[vect_Ugf_z,repmat(vect_Ugf_z(:,end),[1 size(vect_Ugf_z_twin,2)])+vect_Ugf_z_twin];
                Dfft_gruppo_vect=[Dfft_gruppo_vect,repmat(Dfft_gruppo_vect(:,end),[1 size(Dfft_gruppo_vect_twin,2)])+Dfft_gruppo_vect_twin];
                I_gcurve_vect=[I_gcurve_vect,I_gcurve_vect_twin];
                num_incr_TUNNEL=[num_incr_TUNNEL+num_incr_TUNNEL_twin];
                
                ys_vect=[ys_vect,ys_vect_twin+ys_vect(end)-ys_vect_twin(1)];
                
                ys_tun1st=zeros(size(ys_vect));
                ys_tun1st(1,1:num_incr_tunadv)=ys_vect(1,1:num_incr_tunadv);
                ys_tun1st(1,1+num_incr_tunadv:end)=ys_vect(num_incr_tunadv);
                ys_tun2nd=zeros(size(ys_vect));
                ys_tun2nd(1,1+num_incr_tunadv:end)=ys_vect_twin;

                
                %I need to project the displacements from the twin refer system to
                %the first tunnel global reference system
                dxbg_twin_proj=+dxbg_twin*cos(theta-theta_twin)-dybg_twin*sin(theta-theta_twin);
                dybg_twin_proj=+dxbg_twin*sin(theta-theta_twin)+dybg_twin*cos(theta-theta_twin);
                dzbg_twin_proj=dzbg_twin;
              
                dxbg=dxbg+dxbg_twin_proj;
                dybg=dybg+dybg_twin_proj;
                dzbg=dzbg+dzbg_twin_proj;

                
                
                
                
                dxbl=dxbl+dxbl_twin;
                dybl=dybl+dybl_twin;
                dzbl=dzbl+dzbl_twin;
                
                
                
            end
            %*************


else
    

Vltp_vect=0:Vltp/num_incr_TUNNEL:Vltp;
epsilon_vect=Vltp_vect/100/2;

for ii=1:length(Vltp_vect)
    if switch_soil==1
        [ Uffx(:,:,ii),Uffz(:,:,ii) ]=u_LON(z_global_foot,x_global_foot,ht,Rt,epsilon_vect(ii),0,0.4999999);
    elseif switch_soil==101
        [ Uffx(:,:,ii),Uffz(:,:,ii) ]=u_DECK(z_global_foot,x_global_foot,B_X_foot,Rt,switch_soil);
        Uffx(:,:,ii)=Uffx(:,:,ii).*(Vltp_vect(ii)/max(Vltp_vect));
        Uffz(:,:,ii)=Uffz(:,:,ii).*(Vltp_vect(ii)/max(Vltp_vect));
        I_gcurve_vect = B_X_foot*ones(size(Vltp_vect));
    elseif switch_soil==100 
        [ Uffx(:,:,ii),Uffz(:,:,ii) ]=u_DECK(z_global_foot,x_global_foot,B_X_foot,Rt,switch_soil);
        Uffx(:,:,ii)=Uffx(:,:,ii).*(Vltp_vect(ii)/max(Vltp_vect));
        Uffz(:,:,ii)=Uffz(:,:,ii).*(Vltp_vect(ii)/max(Vltp_vect));
        I_gcurve_vect = B_X_foot*ones(size(Vltp_vect));
    elseif switch_soil==10
        [ Uffx(:,:,ii),Uffz(:,:,ii) ]=u_stgauss             (z_global_foot,x_global_foot,ht,Rt,epsilon_vect(ii));
        I_gcurve_vect = 0.5*ht*ones(size(Vltp_vect));
    elseif switch_soil==11
        [ Uffx(:,:,ii),Uffz(:,:,ii) ]=u_stgauss_vecttuncentr(z_global_foot,x_global_foot,ht,Rt,epsilon_vect(ii));
        I_gcurve_vect = 0.5*ht*ones(size(Vltp_vect));
    elseif switch_soil==2
        load out_u_reg_coeff_CD13_ID90.mat
        load out_u_reg_coeffx_CD13_ID90.mat
        id=0.9;
        [ Uffx(:,:,ii),Uffz(:,:,ii) ]=u_SE(z_global_foot,x_global_foot,ht,Rt,epsilon_vect(ii),id,reg_ca,reg_cb,reg_c1,reg_c2,reg_c3,reg_c4,reg_c5,reg_c6,reg_cax,reg_cbx,reg_c1x,reg_c2x,reg_c3x,reg_c4x,reg_c5x,reg_c6x)  ; 
    elseif switch_soil==224
        load out_u_reg_coeff_CD24_ID90_conf2015.mat
        load out_u_reg_coeffx_CD24_ID90_conf2015.mat
        id=0.9;
        I_gcurve_vect = 0.5*ht*ones(size(Vltp_vect));
        [ Uffx(:,:,ii),Uffz(:,:,ii) ]=u_SE(z_global_foot,x_global_foot,ht,Rt,epsilon_vect(ii),id,reg_ca,reg_cb,reg_c1,reg_c2,reg_c3,reg_c4,reg_c5,reg_c6,reg_cax,reg_cbx,reg_c1x,reg_c2x,reg_c3x,reg_c4x,reg_c5x,reg_c6x)  ; 

        
    elseif switch_soil==3000
        [ Uffx(:,:,ii),Uffz(:,:,ii) ]=[ dxbg,-dzbl ] ;
        
    %Centrifuge test data in sand
    elseif switch_soil==3 %Farrell test
        load WS_int_vltmax_CD13_id09.mat
        N=75;
        Rt=N*0.041;
        ht=N*0.150;
        I_gcurve=(I_mg_cd13_id09*N/1000);
        I_gcurve_vect = interp1(Vlt_cd13_id09,I_gcurve,Vltp_vect);
        if ii==1
            Uffz(:,:,ii)=0*x_global_foot;
            Uffx(:,:,ii)=0*x_global_foot;
        else
            Uffz(:,:,ii)=N/1000* griddata(cell2mat(X_int_cd13_id09),cell2mat(Vlt_int_cd13_id09),cell2mat(DY_int_cd13_id09),x_global_foot*1000/N,Vltp_vect(ii)*ones(size(x_global_foot*1000/N)),'linear');
            Uffx(:,:,ii)=N/1000* griddata(cell2mat(X_int_cd13_id09),cell2mat(Vlt_int_cd13_id09),cell2mat(DX_int_cd13_id09),x_global_foot*1000/N,Vltp_vect(ii)*ones(size(x_global_foot*1000/N)),'linear');
        end

elseif switch_soil==30068 %Franza test
        load WS_int_vltmax_CD13_id03.mat
        N=68;
        Rt=N*0.09/2;
        ht=N*0.165;
        I_gcurve=(I_mg_cd13_id03*N/1000);
        I_gcurve_vect = interp1(Vlt_cd13_id03,I_gcurve,Vltp_vect);
        if ii==1
            Uffz(:,:,ii)=0*x_global_foot;
            Uffx(:,:,ii)=0*x_global_foot;
        else
            Uffz(:,:,ii)=N/1000* griddata(cell2mat(X_int_cd13_id03),cell2mat(Vlt_int_cd13_id03),cell2mat(DY_int_cd13_id03),x_global_foot*1000/N,Vltp_vect(ii)*ones(size(x_global_foot*1000/N)),'linear');
            Uffx(:,:,ii)=N/1000* griddata(cell2mat(X_int_cd13_id03),cell2mat(Vlt_int_cd13_id03),cell2mat(DX_int_cd13_id03),x_global_foot*1000/N,Vltp_vect(ii)*ones(size(x_global_foot*1000/N)),'linear');
        end
        
        
    elseif switch_soil==32030
        load WS_int_vltmax_CD20_id03c.mat
        N=80;
        Rt=N*0.09/2;
        ht=N*0.225;
        I_gcurve=(I_mg_cd20_id03c*N/1000);
        I_gcurve_vect = interp1(Vlt_cd20_id03c(Vlt_cd20_id03c~=0),I_gcurve(I_gcurve~=0),Vltp_vect);
        if ii==1
            Uffz(:,:,ii)=0*x_global_foot;
            Uffx(:,:,ii)=0*x_global_foot;
        else
            Uffz(:,:,ii)=N/1000* griddata(cell2mat(X_int_cd20_id03c),cell2mat(Vlt_int_cd20_id03c),cell2mat(DY_int_cd20_id03c),x_global_foot*1000/N,Vltp_vect(ii)*ones(size(x_global_foot*1000/N)),'linear');
            Uffx(:,:,ii)=0*x_global_foot;
        end
        
    elseif switch_soil==32050
        load WS_int_vltmax_CD20_id05.mat
        N=80;
        Rt=N*0.09/2;
        ht=N*0.225;
        I_gcurve=(I_mg_cd20_id05*N/1000);
        I_gcurve_vect = interp1(Vlt_cd20_id05(Vlt_cd20_id05~=0),I_gcurve(I_gcurve~=0),Vltp_vect);
        if ii==1
            Uffz(:,:,ii)=0*x_global_foot;
            Uffx(:,:,ii)=0*x_global_foot;
        else
            Uffz(:,:,ii)=N/1000* griddata(cell2mat(X_int_cd20_id05),cell2mat(Vlt_int_cd20_id05),cell2mat(DY_int_cd20_id05),x_global_foot*1000/N,Vltp_vect(ii)*ones(size(x_global_foot*1000/N)),'linear');
            Uffx(:,:,ii)=0*x_global_foot;
        end   

    elseif switch_soil==32090
        load WS_int_vltmax_CD20_id09.mat
        N=80;
        Rt=N*0.09/2;
        ht=N*0.225;
        I_gcurve=(I_mg_cd20_id09*N/1000);
        I_gcurve_vect = interp1(Vlt_cd20_id09(Vlt_cd20_id09~=0),I_gcurve(I_gcurve~=0),Vltp_vect);
        if ii==1
            Uffz(:,:,ii)=0*x_global_foot;
            Uffx(:,:,ii)=0*x_global_foot;
        else
            Uffz(:,:,ii)=N/1000* griddata(cell2mat(X_int_cd20_id09),cell2mat(Vlt_int_cd20_id09),cell2mat(DY_int_cd20_id09),x_global_foot*1000/N,Vltp_vect(ii)*ones(size(x_global_foot*1000/N)),'linear');
            Uffx(:,:,ii)=0*x_global_foot;
        end   
        
    elseif switch_soil==36330
        load WS_int_vltmax_CD63_id03.mat
        N=80;
        Rt=N*0.04/2;
        ht=N*0.27;
        I_gcurve=(I_mg_cd63_id03*N/1000);
        I_gcurve_vect = interp1(Vlt_cd63_id03(Vlt_cd63_id03~=0),I_gcurve(I_gcurve~=0),Vltp_vect);
        if ii==1
            Uffz(:,:,ii)=0*x_global_foot;
            Uffx(:,:,ii)=0*x_global_foot;
        else
            Uffz(:,:,ii)=N/1000* griddata(cell2mat(X_int_cd63_id03),cell2mat(Vlt_int_cd63_id03),cell2mat(DY_int_cd63_id03),x_global_foot*1000/N,Vltp_vect(ii)*ones(size(x_global_foot*1000/N)),'linear');
            Uffx(:,:,ii)=0*x_global_foot;
        end   
        
    elseif switch_soil==36350
        load WS_int_vltmax_CD63_id05.mat
        N=80;
        Rt=N*0.04/2;
        ht=N*0.27;
        I_gcurve=(I_mg_cd63_id05*N/1000);
        I_gcurve_vect = interp1(Vlt_cd63_id05(Vlt_cd63_id05~=0),I_gcurve(I_gcurve~=0),Vltp_vect);
        if ii==1
            Uffz(:,:,ii)=0*x_global_foot;
            Uffx(:,:,ii)=0*x_global_foot;
        else
            Uffz(:,:,ii)=N/1000* griddata(cell2mat(X_int_cd63_id05),cell2mat(Vlt_int_cd63_id05),cell2mat(DY_int_cd63_id05),x_global_foot*1000/N,Vltp_vect(ii)*ones(size(x_global_foot*1000/N)),'linear');
            Uffx(:,:,ii)=0*x_global_foot;
        end   
        
     elseif switch_soil==36390
        load WS_int_vltmax_CD63_id09.mat
        N=80;
        Rt=N*0.04/2;
        ht=N*0.27;
        I_gcurve=(I_mg_cd63_id09*N/1000);
        I_gcurve_vect = interp1(Vlt_cd63_id09(Vlt_cd63_id09~=0),I_gcurve(I_gcurve~=0),Vltp_vect);
        if ii==1
            Uffz(:,:,ii)=0*x_global_foot;
            Uffx(:,:,ii)=0*x_global_foot;
        else
            Uffz(:,:,ii)=N/1000* griddata(cell2mat(X_int_cd63_id09),cell2mat(Vlt_int_cd63_id09),cell2mat(DY_int_cd63_id09),x_global_foot*1000/N,Vltp_vect(ii)*ones(size(x_global_foot*1000/N)),'linear');
            Uffx(:,:,ii)=0*x_global_foot;
        end   
        
    end
        
    Uffx_vect(:,:,ii)=coeffux*Uffx(:,:,ii);
    Uffz_vect(:,:,ii)=coeffuz*Uffz(:,:,ii);
    Uffy_vect(:,:,ii)=zeros(nnodes_foot,nfoot);

   
    Uffz(isnan(Uffz)) = 0 ;
	Uffx(isnan(Uffx)) = 0 ;
	Uffx_vect(isnan(Uffx_vect)) = 0 ;
	Uffz_vect(isnan(Uffz_vect)) = 0 ;
    

    
end




%Uff VECTOR AT NODES OF FINITE ELEMENTS
Dfft_gruppo_vect=zeros(nfoot*nnodes_foot*6,length(Vltp_vect));
for ii=1:length(Vltp_vect)
    for i=1:nnodes_foot;
        for ip=1:nfoot;
            Dfft_gruppo_vect((ip-1)*(nnodes_foot*6)+(i-1)*6+1,ii)=Uffx_vect(i,ip,ii);
            Dfft_gruppo_vect((ip-1)*(nnodes_foot*6)+(i-1)*6+2,ii)=Uffy_vect(i,ip,ii);
            Dfft_gruppo_vect((ip-1)*(nnodes_foot*6)+(i-1)*6+3,ii)=Uffz_vect(i,ip,ii);
        end
    end
end

vect_Ugf_z=Dfft_gruppo_vect(3:6:end,:);
vect_Ugf_x=Dfft_gruppo_vect(1:6:end,:);


end



%% 
% =========================================================================
% ELASTO SOLUTION (KLAR ET AL 2005)
% =========================================================================


%STIFFNESS AND FLEXIBILITY MATRICES OF THE SOIL (Only translational dofs)
lamdasts_rid=FLEX-diag(diag(FLEX));
lamdastd_rid=diag(diag(FLEX));
Kst_rid=inv(FLEX-lamdasts_rid);

% Zeros are added at the rotational dofs
Kst=zeros(6*nfoot*(nelementi_foot+1));
lamdasts=zeros(6*nfoot*(nelementi_foot+1));
lamdastd=zeros(6*nfoot*(nelementi_foot+1));
for i=1:nfoot*(nelementi_foot+1)
    for j=1:nfoot*(nelementi_foot+1)
        Kst(1+(i-1)*6:3+(i-1)*6,1+(j-1)*6:3+(j-1)*6)=...
            Kst_rid(1+(i-1)*3:3+(i-1)*3,1+(j-1)*3:3+(j-1)*3);
        lamdasts(1+(i-1)*6:3+(i-1)*6,1+(j-1)*6:3+(j-1)*6)=...
            lamdasts_rid(1+(i-1)*3:3+(i-1)*3,1+(j-1)*3:3+(j-1)*3);
        lamdastd(1+(i-1)*6:3+(i-1)*6,1+(j-1)*6:3+(j-1)*6)=...
            lamdastd_rid(1+(i-1)*3:3+(i-1)*3,1+(j-1)*3:3+(j-1)*3);
    end
end


%TOTAL STIFFNESS MATRIX
KK_el=KKfoot+Kst+KKextra+Kst*lamdasts*KKfoot;

KK_fr_i_el=(KKfoot+KKframe_fixed_foot)+Kst+KKextra+Kst*lamdasts*(KKfoot+KKframe_fixed_foot);
KK_fr_h_el=(KKfoot+KKframe_hinged_foot)+Kst+KKextra+Kst*lamdasts*(KKfoot+KKframe_hinged_foot);




for jj=1:num_incr_TUNNEL
    
    %GF INPUT
    UCAT=Dfft_gruppo_vect(:,jj+1); %Tunnelling-induced greenfield movements
        
    %FORCE VECTOR
    %Tunnelling-induced forces
    F_el=Kst*UCAT;
    
    %External forces (concentrated at the centre of the footings)
    P_el=zeros(size(F_el));
    
    P_el(1+(nelementi_foot/2*6):6*nnodes_foot:end,1)=Fx_foot_centr;
    P_el(3+(nelementi_foot/2*6):6*nnodes_foot:end,1)=Fz_foot_centr;
    P_el(5+(nelementi_foot/2*6):6*nnodes_foot:end,1)=My_foot_centr;


    
   
    %External distributed loads along the footings axis
    q_ind=ones(size(P_el));
    q_ind(1:6*nnodes_foot:end,1)=0.5;
    q_ind(3:6*nnodes_foot:end,1)=0.5;
    q_ind(1+(nelementi_foot*6):6*nnodes_foot:end,1)=0.5;
    q_ind(3+(nelementi_foot*6):6*nnodes_foot:end,1)=0.5;
    if length(qx_foot)>1
        error('qx and qz should have a unique value')
    end
    P_el(1:6:end,1)                     =q_ind(1:6:end,1).*qx_foot(1)*h_el_foot+P_el(1:6:end,1);
    P_el(3:6:end,1)                     =q_ind(3:6:end,1).*qz_foot(1)*h_el_foot+P_el(3:6:end,1);

    
    %External distributed loads along the beam of the frame

	if length(X_foot_centr)>1
	Fwx_foot_centr=Reactions_qbeam_fixhead_foot(1:3:end);
    Fwy_foot_centr=Reactions_qbeam_fixhead_foot(2:3:end);
    Mwy_foot_centr=Reactions_qbeam_fixhead_foot(3:3:end);

    W_el=zeros(size(F_el));
	W_el(1+(nelementi_foot/2*6):6*nnodes_foot:end,1)=Fwx_foot_centr;
    W_el(3+(nelementi_foot/2*6):6*nnodes_foot:end,1)=Fwy_foot_centr;
    W_el(5+(nelementi_foot/2*6):6*nnodes_foot:end,1)=Mwy_foot_centr;
    
    elseif length(X_foot_centr)==1 && max(spanx_frame)>0&& max(spanz_frame)>0
    Fwx_foot_centr=Reactions_qbeam_fixhead_foot(1:3:end);
    Fwy_foot_centr=Reactions_qbeam_fixhead_foot(2:3:end);
    Mwy_foot_centr=Reactions_qbeam_fixhead_foot(3:3:end);

    W_el=zeros(size(F_el));
	W_el(1:6*nelementi_bayraft:end,1)=Fwx_foot_centr;
    W_el(3:6*nelementi_bayraft:end,1)=Fwy_foot_centr;
    W_el(5:6*nelementi_bayraft:end,1)=Mwy_foot_centr;  
    else
    W_el=zeros(size(F_el));

    end
    
    P_el=P_el+W_el;   

    
    
    

    

    
    
    
    %SOLUTION DISPLACEMENT VECTOR
    U_node_free_el_T(:,jj)=(KK_el)      \(Kst*lamdasts*P_el+P_el+F_el);
    U_node_fr_i_el_T(:,jj)=(KK_fr_i_el) \(Kst*lamdasts*P_el+P_el+F_el);
    U_node_fr_h_el_T(:,jj)=(KK_fr_h_el) \(Kst*lamdasts*P_el+P_el+F_el);
    
    U_node_free_el_P(:,jj)=(KK_el)      \(Kst*lamdasts*P_el+P_el);
    U_node_fr_i_el_P(:,jj)=(KK_fr_i_el) \(Kst*lamdasts*P_el+P_el);
    U_node_fr_h_el_P(:,jj)=(KK_fr_h_el) \(Kst*lamdasts*P_el+P_el);
    
    F_soil_fr_i_el_T(:,jj)=(KKfoot+KKframe_fixed_foot)*U_node_fr_i_el_T(:,jj)-P_el;
    F_soil_fr_i_el_P(:,jj)=(KKfoot+KKframe_fixed_foot)*U_node_fr_i_el_P(:,jj)-P_el;
    F_soil_fr_i_el_T_z(:,jj)=F_soil_fr_i_el_T(3:6:end,jj);
    F_soil_fr_i_el_T_x(:,jj)=F_soil_fr_i_el_T(1:6:end,jj);

    
    
    F_soil_fr_i_el_deltaT(:,jj)=F_soil_fr_i_el_T(:,jj)-F_soil_fr_i_el_P(:,jj);
	F_soil_fr_i_el_deltaT_z(:,jj)=F_soil_fr_i_el_deltaT(3:6:end,jj);
    F_soil_fr_i_el_deltaT_x(:,jj)=F_soil_fr_i_el_deltaT(1:6:end,jj);
    
    
end


vect_Utotf_fr_i_el_P=U_node_fr_i_el_P;
vect_Utotf_fr_i_el_T=U_node_fr_i_el_T;
vect_Utotp_fr_i_el_T=zeros(size(U_node_fr_i_el_P));
vect_Utotp_fr_i_el_P=zeros(size(U_node_fr_i_el_P));

for i=1:nnodes_foot;
    Uxpg_fr_i_el_P(i,1:1:nfoot)  =vect_Utotf_fr_i_el_P(1+(i-1)*6:nnodes_foot*6:end,end);
    Uypg_fr_i_el_P(i,1:1:nfoot)  =vect_Utotf_fr_i_el_P(2+(i-1)*6:nnodes_foot*6:end,end);
    Uzpg_fr_i_el_P(i,1:1:nfoot)  =vect_Utotf_fr_i_el_P(3+(i-1)*6:nnodes_foot*6:end,end);
    Ur2pg_fr_i_el_P(i,1:1:nfoot) =vect_Utotf_fr_i_el_P(5+(i-1)*6:nnodes_foot*6:end,end);
    
    Uxpg_fr_i_el_T(i,1:1:nfoot)  =vect_Utotf_fr_i_el_T(1+(i-1)*6:nnodes_foot*6:end,end);
    Uypg_fr_i_el_T(i,1:1:nfoot)  =vect_Utotf_fr_i_el_T(2+(i-1)*6:nnodes_foot*6:end,end);
    Uzpg_fr_i_el_T(i,1:1:nfoot)  =vect_Utotf_fr_i_el_T(3+(i-1)*6:nnodes_foot*6:end,end);
    Ur2pg_fr_i_el_T(i,1:1:nfoot) =vect_Utotf_fr_i_el_T(5+(i-1)*6:nnodes_foot*6:end,end);
end

Uxpg_fr_i_el_deltaT  =Uxpg_fr_i_el_T-Uxpg_fr_i_el_P;
Uypg_fr_i_el_deltaT  =Uypg_fr_i_el_T-Uypg_fr_i_el_P;
Uzpg_fr_i_el_deltaT  =Uzpg_fr_i_el_T-Uzpg_fr_i_el_P;
Ur2pg_fr_i_el_deltaT =Ur2pg_fr_i_el_T-Ur2pg_fr_i_el_P;

vect_Utotp_fr_i_el_P_z=vect_Utotp_fr_i_el_P(3:6:end,:);
vect_Utotp_fr_i_el_P_x=vect_Utotp_fr_i_el_P(1:6:end,:);
vect_Utotf_fr_i_el_P_z=vect_Utotf_fr_i_el_P(3:6:end,:);
vect_Utotf_fr_i_el_P_x=vect_Utotf_fr_i_el_P(1:6:end,:);
vect_Utotp_fr_i_el_T_z=vect_Utotp_fr_i_el_T(3:6:end,:);
vect_Utotp_fr_i_el_T_x=vect_Utotp_fr_i_el_T(1:6:end,:);
vect_Utotf_fr_i_el_T_z=vect_Utotf_fr_i_el_T(3:6:end,:);
vect_Utotf_fr_i_el_T_x=vect_Utotf_fr_i_el_T(1:6:end,:);




% vect_Utotf_fr_i_el_deltaT_z=vect_Utotf_fr_i_el_T_z-vect_Utotf_fr_i_el_P_z(:,end);
% vect_Utotf_fr_i_el_deltaT_x=vect_Utotf_fr_i_el_T_x-vect_Utotf_fr_i_el_P_x(:,end);

vect_Utotp_fr_i_el_P_r2=vect_Utotp_fr_i_el_P(5:6:end,:);
vect_Utotf_fr_i_el_P_r2=vect_Utotf_fr_i_el_P(5:6:end,:);
vect_Utotp_fr_i_el_T_r2=vect_Utotp_fr_i_el_T(5:6:end,:);
vect_Utotf_fr_i_el_T_r2=vect_Utotf_fr_i_el_T(5:6:end,:);

vect_Utotf_fr_i_el_deltaT_z = vect_Utotf_fr_i_el_T_z - repmat(vect_Utotf_fr_i_el_P_z(:,end), [1, size(vect_Utotf_fr_i_el_T_z,2)]);
vect_Utotf_fr_i_el_deltaT_x = vect_Utotf_fr_i_el_T_x - repmat(vect_Utotf_fr_i_el_P_x(:,end), [1, size(vect_Utotf_fr_i_el_T_x,2)]);
vect_Utotf_fr_i_el_deltaT_r2 = vect_Utotf_fr_i_el_T_r2 - repmat(vect_Utotf_fr_i_el_P_r2(:,end), [1, size(vect_Utotf_fr_i_el_T_r2,2)]);

vect_Utotf_fr_i_el_deltaT = vect_Utotf_fr_i_el_T - repmat(vect_Utotf_fr_i_el_P(:,end), [1, size(vect_Utotf_fr_i_el_T,2)]);


 %% 
% =========================================================================
% DISPLACEMENT SOLUTION according to FRANZA ET AL 2017
% =========================================================================





%When I contrain rotation I obtain horizontal force
Aext=zeros(size(KKsoil));
aext=zeros(6,6);
aext(1,5)=1*d_na;%assumption small displacement tan(phi)=phi thus ux=d/2*phi
for ij=1:size(Aext,1)/6
Aext(1+(ij-1)*6:6+(ij-1)*6,1+(ij-1)*6:6+(ij-1)*6)=aext;
end
%When I contrain rotation I obtain bending moment from horizontal force
Bext=zeros(size(KKsoil));
bext=zeros(6,6);
bext(5,1)=1*d_na;
for ij=1:size(Bext,1)/6
Bext(1+(ij-1)*6:6+(ij-1)*6,1+(ij-1)*6:6+(ij-1)*6)=bext;
end

Kcoupl=(+KKsoil*Aext)'+KKsoil*Aext+Bext*KKsoil*Aext;


  Kcouplv2=  KKsoil;
  for ii=1:size(Kcouplv2,1)/6
  Kcouplv2(:,5+(ii-1)*6)= d_na* Kcouplv2(:,1+(ii-1)*6);
  end
  
  for ii=1:size(Kcouplv2,1)/6
  Kcouplv2(5+(ii-1)*6,:)= d_na* Kcouplv2(1+(ii-1)*6,:);
  end
  
  Kcouplv2=Kcouplv2-KKsoil;
  
BBcheck=Kcouplv2-Kcoupl;
BB=Kcoupl-Kcoupl';


%EXTERNAL FORCES DUE TO TUNNELLING
for jj=2:num_incr_TUNNEL+1
ii=jj-1;
Fpg=(KKsoil)*Dfft_gruppo_vect(:,jj); %vettore forze esterne
U_node_fr_i(:,ii)=(KKpg_fr_i)\(P_el+Fpg);
U_node_fr_i_P(:,ii)=(KKpg_fr_i)\(P_el+0);
U_node_fr_i_deltaT(:,ii)=(KKpg_fr_i)\(Fpg);
U_node_fr_h(:,ii)=(KKpg_fr_h)\(P_el+Fpg); 
U_node_free(:,ii)=(KKpg)\(P_el+Fpg); 
Fpg_GL=Fpg; %vettore forze esterne at the ground level
Fpg_GL_z=Fpg_GL(3:6:end);Fpg_GL_x=Fpg_GL(1:6:end);


Fpg(5:6:end)=Fpg(1:6:end)*(1*d_na);
U_node_fr_i_off(:,ii)=(KKpg_fr_i+Kcoupl)\(P_el+Fpg);
U_node_fr_i_off_deltaT(:,ii)=(KKpg_fr_i+Kcoupl)\(Fpg);
U_node_fr_i_off_P(:,ii)=(KKpg_fr_i+Kcoupl)\(P_el);
end



for jj=2:num_incr_TUNNEL+1
ii=jj-1;
for i=1:nnodes_foot;
	for ip=1:nfoot;
% 	Fxpg(i,ip)=Fpg((ip-1)*(nnodes_foot*6)+((i-1)*6)+1);
%     Fypg(i,ip)=Fpg((ip-1)*(nnodes_foot*6)+((i-1)*6)+2);
%     Fzpg(i,ip)=Fpg((ip-1)*(nnodes_foot*6)+((i-1)*6)+3);
       
%     Uxpg_fr_i(i,ip,ii)=U_node_fr_i((ip-1)*(nnodes_foot*6)+((i-1)*6)+1,ii);
%     Uypg_fr_i(i,ip,ii)=U_node_fr_i((ip-1)*(nnodes_foot*6)+((i-1)*6)+2,ii);
%     Uzpg_fr_i(i,ip,ii)=U_node_fr_i((ip-1)*(nnodes_foot*6)+((i-1)*6)+3,ii);
%     Ur2pg_fr_i(i,ip,ii)=U_node_fr_i((ip-1)*(nnodes_foot*6)+((i-1)*6)+5,ii);
% 
%     Uxpg_fr_i_deltaT(i,ip,ii)=U_node_fr_i_deltaT((ip-1)*(nnodes_foot*6)+((i-1)*6)+1,ii);
%     Uypg_fr_i_deltaT(i,ip,ii)=U_node_fr_i_deltaT((ip-1)*(nnodes_foot*6)+((i-1)*6)+2,ii);
%     Uzpg_fr_i_deltaT(i,ip,ii)=U_node_fr_i_deltaT((ip-1)*(nnodes_foot*6)+((i-1)*6)+3,ii);
%     Ur2pg_fr_i_deltaT(i,ip,ii)=U_node_fr_i_deltaT((ip-1)*(nnodes_foot*6)+((i-1)*6)+5,ii);    
%     
%     
%     Uxpg_fr_i_off(i,ip,ii)=U_node_fr_i_off((ip-1)*(nnodes_foot*6)+((i-1)*6)+1,ii);
%     Uypg_fr_i_off(i,ip,ii)=U_node_fr_i_off((ip-1)*(nnodes_foot*6)+((i-1)*6)+2,ii);
%     Uzpg_fr_i_off(i,ip,ii)=U_node_fr_i_off((ip-1)*(nnodes_foot*6)+((i-1)*6)+3,ii);
%     Ur2pg_fr_i_off(i,ip,ii)=U_node_fr_i_off((ip-1)*(nnodes_foot*6)+((i-1)*6)+5,ii);
% 
%     Uxpg_fr_i_off_deltaT(i,ip,ii)=U_node_fr_i_off_deltaT((ip-1)*(nnodes_foot*6)+((i-1)*6)+1,ii);
%     Uypg_fr_i_off_deltaT(i,ip,ii)=U_node_fr_i_off_deltaT((ip-1)*(nnodes_foot*6)+((i-1)*6)+2,ii);
%     Uzpg_fr_i_off_deltaT(i,ip,ii)=U_node_fr_i_off_deltaT((ip-1)*(nnodes_foot*6)+((i-1)*6)+3,ii);
%     Ur2pg_fr_i_off_deltaT(i,ip,ii)=U_node_fr_i_off_deltaT((ip-1)*(nnodes_foot*6)+((i-1)*6)+5,ii);    
%     
%     Uxpg_fr_h(i,ip,ii)=U_node_fr_h((ip-1)*(nnodes_foot*6)+((i-1)*6)+1,ii);
%     Uypg_fr_h(i,ip,ii)=U_node_fr_h((ip-1)*(nnodes_foot*6)+((i-1)*6)+2,ii);
%     Uzpg_fr_h(i,ip,ii)=U_node_fr_h((ip-1)*(nnodes_foot*6)+((i-1)*6)+3,ii);   
%     Ur2pg_fr_h(i,ip,ii)=U_node_fr_h((ip-1)*(nnodes_foot*6)+((i-1)*6)+5,ii);
%     
%     Uxpg_free(i,ip,ii)=U_node_free((ip-1)*(nnodes_foot*6)+((i-1)*6)+1,ii);
%     Uypg_free(i,ip,ii)=U_node_free((ip-1)*(nnodes_foot*6)+((i-1)*6)+2,ii);
%     Uzpg_free(i,ip,ii)=U_node_free((ip-1)*(nnodes_foot*6)+((i-1)*6)+3,ii);
%     Ur2pg_free(i,ip,ii)=U_node_free((ip-1)*(nnodes_foot*6)+((i-1)*6)+5,ii);

    Uxpg_fr_i(:,ii)=U_node_fr_i(1:6:end,ii);
    Uypg_fr_i(:,ii)=U_node_fr_i(2:6:end,ii);
    Uzpg_fr_i(:,ii)=U_node_fr_i(3:6:end,ii);
    Ur2pg_fr_i(:,ii)=U_node_fr_i(5:6:end,ii);

    Uxpg_fr_i_deltaT(:,ii)=U_node_fr_i_deltaT(1:6:end,ii);
    Uypg_fr_i_deltaT(:,ii)=U_node_fr_i_deltaT(2:6:end,ii);
    Uzpg_fr_i_deltaT(:,ii)=U_node_fr_i_deltaT(3:6:end,ii);
    Ur2pg_fr_i_deltaT(:,ii)=U_node_fr_i_deltaT(5:6:end,ii);    
    
    
    Uxpg_fr_i_off(:,ii)=U_node_fr_i_off(1:6:end,ii);
    Uypg_fr_i_off(:,ii)=U_node_fr_i_off(2:6:end,ii);
    Uzpg_fr_i_off(:,ii)=U_node_fr_i_off(3:6:end,ii);
    Ur2pg_fr_i_off(:,ii)=U_node_fr_i_off(5:6:end,ii);

    Uxpg_fr_i_off_deltaT(:,ii)=U_node_fr_i_off_deltaT(1:6:end,ii);
    Uypg_fr_i_off_deltaT(:,ii)=U_node_fr_i_off_deltaT(2:6:end,ii);
    Uzpg_fr_i_off_deltaT(:,ii)=U_node_fr_i_off_deltaT(3:6:end,ii);
    Ur2pg_fr_i_off_deltaT(:,ii)=U_node_fr_i_off_deltaT(5:6:end,ii);    
    
    Uxpg_fr_h(:,ii)=U_node_fr_h(1:6:end,ii);
    Uypg_fr_h(:,ii)=U_node_fr_h(2:6:end,ii);
    Uzpg_fr_h(:,ii)=U_node_fr_h(3:6:end,ii);   
    Ur2pg_fr_h(:,ii)=U_node_fr_h(5:6:end,ii);
    
    Uxpg_free(:,ii)=U_node_free(1:6:end,ii);
    Uypg_free(:,ii)=U_node_free(2:6:end,ii);
    Uzpg_free(:,ii)=U_node_free(3:6:end,ii);
    Ur2pg_free(:,ii)=U_node_free(5:6:end,ii);


    end
end
end


%%
% =========================================================================
% ELASTO-PLASTIC SOLUTION (KLAR ET AL 2007)
% =========================================================================

if switch_ep==0
vect_Utotf_fr_i_ep_P=vect_Utotf_fr_i_el_P;
vect_Utotf_fr_i_ep_T=vect_Utotf_fr_i_el_T;
vect_Utotf_fr_i_ep_deltaT_z=vect_Utotf_fr_i_el_deltaT_z;
vect_Utotf_fr_i_ep_deltaT_x=vect_Utotf_fr_i_el_deltaT_x;
vect_Utotf_fr_i_ep_deltaT_r2=vect_Utotf_fr_i_el_deltaT_r2;   

vect_Utotf_fr_i_ep_deltaT=vect_Utotf_fr_i_el_deltaT;
end



if switch_ep==1

% ------------------------------------------------------------------------- 
% FREE FOOTINGS
% -------------------------------------------------------------------------

SS=KKfoot+KKextra;

script_ep_solution

vect_Utotf_free_ep_P=vect_Utotf_P;
vect_Utotf_free_ep_T=vect_Utotf_T;
vect_Utotp_free_ep_T=cumsum(vect_Upinc_T,2);
vect_Utotp_free_ep_P=cumsum(vect_Upinc_P,2);

F_soil_free_ep=SS*vect_Utotf_T(:,end)-P_el;
F_soil_free_ep_z=F_soil_free_ep(3:6:end);
F_soil_free_ep_x=F_soil_free_ep(1:6:end);

for i=1:nnodes_foot;
    Uxpg_free_ep_P(i,1:1:nfoot)  =vect_Utotf_P(1+(i-1)*6:nnodes_foot*6:end,end);
    Uypg_free_ep_P(i,1:1:nfoot)  =vect_Utotf_P(2+(i-1)*6:nnodes_foot*6:end,end);
    Uzpg_free_ep_P(i,1:1:nfoot)  =vect_Utotf_P(3+(i-1)*6:nnodes_foot*6:end,end);
    Ur2pg_free_ep_P(i,1:1:nfoot) =vect_Utotf_P(5+(i-1)*6:nnodes_foot*6:end,end);
    
    Uxpg_free_ep_T(i,1:1:nfoot)  =vect_Utotf_T(1+(i-1)*6:nnodes_foot*6:end,end);
    Uypg_free_ep_T(i,1:1:nfoot)  =vect_Utotf_T(2+(i-1)*6:nnodes_foot*6:end,end);
    Uzpg_free_ep_T(i,1:1:nfoot)  =vect_Utotf_T(3+(i-1)*6:nnodes_foot*6:end,end);
    Ur2pg_free_ep_T(i,1:1:nfoot) =vect_Utotf_T(5+(i-1)*6:nnodes_foot*6:end,end);    
end

Uxpg_free_ep_deltaT  =Uxpg_free_ep_T-Uxpg_free_ep_P;
Uypg_free_ep_deltaT  =Uypg_free_ep_T-Uypg_free_ep_P;
Uzpg_free_ep_deltaT  =Uzpg_free_ep_T-Uzpg_free_ep_P;
Ur2pg_free_ep_deltaT =Ur2pg_free_ep_T-Ur2pg_free_ep_P;

vect_Utotp_free_ep_P_z=vect_Utotp_free_ep_P(3:6:end,:);
vect_Utotp_free_ep_P_x=vect_Utotp_free_ep_P(1:6:end,:);
vect_Utotf_free_ep_P_z=vect_Utotf_free_ep_P(3:6:end,:);
vect_Utotf_free_ep_P_x=vect_Utotf_free_ep_P(1:6:end,:);
vect_Utotp_free_ep_T_z=vect_Utotp_free_ep_T(3:6:end,:);
vect_Utotp_free_ep_T_x=vect_Utotp_free_ep_T(1:6:end,:);
vect_Utotf_free_ep_T_z=vect_Utotf_free_ep_T(3:6:end,:);
vect_Utotf_free_ep_T_x=vect_Utotf_free_ep_T(1:6:end,:);

% vect_Utotf_free_ep_deltaT_z=vect_Utotf_free_ep_T_z-vect_Utotf_free_ep_P_z(:,end);
% vect_Utotf_free_ep_deltaT_x=vect_Utotf_free_ep_T_x-vect_Utotf_free_ep_P_x(:,end);

vect_Utotf_free_ep_deltaT_z = vect_Utotf_free_ep_T_z - repmat(vect_Utotf_free_ep_P_z(:,end), [1, size(vect_Utotf_free_ep_T_z,2)]);
vect_Utotf_free_ep_deltaT_x = vect_Utotf_free_ep_T_x - repmat(vect_Utotf_free_ep_P_x(:,end), [1, size(vect_Utotf_free_ep_T_x,2)]);


%
% ------------------------------------------------------------------------- 
% FRAME ON FOOTINGS
% -------------------------------------------------------------------------

SS=+KKextra+KKfoot+KKframe_fixed_foot;

script_ep_solution

vect_Utotf_fr_i_ep_P=vect_Utotf_P;
vect_Utotf_fr_i_ep_T=vect_Utotf_T;
vect_Utotp_fr_i_ep_T=cumsum(vect_Upinc_T,2);
vect_Utotp_fr_i_ep_P=cumsum(vect_Upinc_P,2);





for ii=1:num_incr_TUNNEL
    F_soil_fr_i_ep_T(:,ii)=SS*vect_Utotf_T(:,ii)-P_el;
end
F_soil_fr_i_ep_T_z=F_soil_fr_i_ep_T(3:6:end,:);
F_soil_fr_i_ep_T_x=F_soil_fr_i_ep_T(1:6:end,:);





sigma_soil_fr_i_el_T_z=F_soil_fr_i_el_T_z./repmat(el_area,[1,num_incr_TUNNEL]);
sigma_soil_fr_i_el_T_x=F_soil_fr_i_el_T_x./repmat(el_area,[1,num_incr_TUNNEL]);
sigma_soil_fr_i_ep_T_z=F_soil_fr_i_ep_T_z./repmat(el_area,[1,num_incr_TUNNEL]);
sigma_soil_fr_i_ep_T_x=F_soil_fr_i_ep_T_x./repmat(el_area,[1,num_incr_TUNNEL]);

sigma_soil_fr_i_el_deltaT_z=F_soil_fr_i_el_deltaT_z./repmat(el_area,[1,num_incr_TUNNEL]);
sigma_soil_fr_i_el_deltaT_x=F_soil_fr_i_el_deltaT_x./repmat(el_area,[1,num_incr_TUNNEL]);




for i=1:nnodes_foot;
    Uxpg_fr_i_ep_P(i,1:1:nfoot)  =vect_Utotf_P(1+(i-1)*6:nnodes_foot*6:end,end);
    Uypg_fr_i_ep_P(i,1:1:nfoot)  =vect_Utotf_P(2+(i-1)*6:nnodes_foot*6:end,end);
    Uzpg_fr_i_ep_P(i,1:1:nfoot)  =vect_Utotf_P(3+(i-1)*6:nnodes_foot*6:end,end);
    Ur2pg_fr_i_ep_P(i,1:1:nfoot) =vect_Utotf_P(5+(i-1)*6:nnodes_foot*6:end,end);
    
    Uxpg_fr_i_ep_T(i,1:1:nfoot)  =vect_Utotf_T(1+(i-1)*6:nnodes_foot*6:end,end);
    Uypg_fr_i_ep_T(i,1:1:nfoot)  =vect_Utotf_T(2+(i-1)*6:nnodes_foot*6:end,end);
    Uzpg_fr_i_ep_T(i,1:1:nfoot)  =vect_Utotf_T(3+(i-1)*6:nnodes_foot*6:end,end);
    Ur2pg_fr_i_ep_T(i,1:1:nfoot) =vect_Utotf_T(5+(i-1)*6:nnodes_foot*6:end,end);    
end

Uxpg_fr_i_ep_deltaT  =Uxpg_fr_i_ep_T-Uxpg_fr_i_ep_P;
Uypg_fr_i_ep_deltaT  =Uypg_fr_i_ep_T-Uypg_fr_i_ep_P;
Uzpg_fr_i_ep_deltaT  =Uzpg_fr_i_ep_T-Uzpg_fr_i_ep_P;
Ur2pg_fr_i_ep_deltaT =Ur2pg_fr_i_ep_T-Ur2pg_fr_i_ep_P;

vect_Utotp_fr_i_ep_P_z=vect_Utotp_fr_i_ep_P(3:6:end,:);
vect_Utotp_fr_i_ep_P_x=vect_Utotp_fr_i_ep_P(1:6:end,:);
vect_Utotf_fr_i_ep_P_z=vect_Utotf_fr_i_ep_P(3:6:end,:);
vect_Utotf_fr_i_ep_P_x=vect_Utotf_fr_i_ep_P(1:6:end,:);
vect_Utotp_fr_i_ep_T_z=vect_Utotp_fr_i_ep_T(3:6:end,:);
vect_Utotp_fr_i_ep_T_x=vect_Utotp_fr_i_ep_T(1:6:end,:);
vect_Utotf_fr_i_ep_T_z=vect_Utotf_fr_i_ep_T(3:6:end,:);
vect_Utotf_fr_i_ep_T_x=vect_Utotf_fr_i_ep_T(1:6:end,:);

vect_Utotp_fr_i_ep_P_r2=vect_Utotp_fr_i_ep_P(5:6:end,:);
vect_Utotf_fr_i_ep_P_r2=vect_Utotf_fr_i_ep_P(5:6:end,:);
vect_Utotp_fr_i_ep_T_r2=vect_Utotp_fr_i_ep_T(5:6:end,:);
vect_Utotf_fr_i_ep_T_r2=vect_Utotf_fr_i_ep_T(5:6:end,:);



% vect_Utotf_fr_i_ep_deltaT_z=vect_Utotf_fr_i_ep_T_z-vect_Utotf_fr_i_ep_P_z(:,end);
% vect_Utotf_fr_i_ep_deltaT_x=vect_Utotf_fr_i_ep_T_x-vect_Utotf_fr_i_ep_P_x(:,end);

vect_Utotf_fr_i_ep_deltaT   = vect_Utotf_fr_i_ep_T   - repmat(vect_Utotf_fr_i_ep_P  (:,end), [1, size(vect_Utotf_fr_i_ep_T  ,2)]);
vect_Utotf_fr_i_ep_deltaT_z = vect_Utotf_fr_i_ep_T_z - repmat(vect_Utotf_fr_i_ep_P_z(:,end), [1, size(vect_Utotf_fr_i_ep_T_z,2)]);
vect_Utotf_fr_i_ep_deltaT_x = vect_Utotf_fr_i_ep_T_x - repmat(vect_Utotf_fr_i_ep_P_x(:,end), [1, size(vect_Utotf_fr_i_ep_T_x,2)]);
vect_Utotf_fr_i_ep_deltaT_r2 = vect_Utotf_fr_i_ep_T_r2 - repmat(vect_Utotf_fr_i_ep_P_r2(:,end), [1, size(vect_Utotf_fr_i_ep_T_r2,2)]);
end
%% FORMATTING VECTORS

Dfft_gruppo_vect(:,1)=[];
vect_Ugf_x(:,1)=[];
vect_Ugf_z(:,1)=[];
Vltp_vect(1)=[];
I_gcurve_vect(1)=[];

if switch_tun_adv_twin==1
Vltp_vect_twin(1)=[];
Vltp_vect=[Vltp_vect,Vltp_vect_twin];

Dfft_gruppo_vect(:,end)=[];
vect_Ugf_x(:,end)=[];
vect_Ugf_z(:,end)=[];

end
%% 

    
    

if length(X_foot_centr)>1
    
    %REACTION FORCES
    
    for i_p=1:nfoot
        for ii=1:size(vect_Utotf_fr_i_ep_T,2)
            M_Utotf_fr_i_ep_deltaT(i_p,ii)=-KKframe_fixed_foot (5+(nelementi_foot/2*6)+6*(nnodes_foot)*(i_p-1),:)*vect_Utotf_fr_i_ep_deltaT(:,ii);
            N_Utotf_fr_i_ep_deltaT(i_p,ii)=-KKframe_fixed_foot (3+(nelementi_foot/2*6)+6*(nnodes_foot)*(i_p-1),:)*vect_Utotf_fr_i_ep_deltaT(:,ii);
            V_Utotf_fr_i_ep_deltaT(i_p,ii)=-KKframe_fixed_foot (1+(nelementi_foot/2*6)+6*(nnodes_foot)*(i_p-1),:)*vect_Utotf_fr_i_ep_deltaT(:,ii);
        end
    end
    
    for i_p=1:nfoot
        M_Utotf_fr_i_ep_P(i_p,1)=-KKframe_fixed_foot (5+(nelementi_foot/2*6)+6*(nnodes_foot)*(i_p-1),:)*vect_Utotf_fr_i_ep_P(:,end);
        N_Utotf_fr_i_ep_P(i_p,1)=-KKframe_fixed_foot (3+(nelementi_foot/2*6)+6*(nnodes_foot)*(i_p-1),:)*vect_Utotf_fr_i_ep_P(:,end);
        V_Utotf_fr_i_ep_P(i_p,1)=-KKframe_fixed_foot (1+(nelementi_foot/2*6)+6*(nnodes_foot)*(i_p-1),:)*vect_Utotf_fr_i_ep_P(:,end);
    end
    
    
    for i_p=1:nfoot
        for ii=1:size(vect_Utotf_fr_i_ep_T,2)
            M_Utotf_fr_i_ep_deltaT(i_p,ii)=-KKframe_fixed_foot (5+(nelementi_foot/2*6)+6*(nnodes_foot)*(i_p-1),:)*vect_Utotf_fr_i_ep_deltaT(:,ii);
            N_Utotf_fr_i_ep_deltaT(i_p,ii)=-KKframe_fixed_foot (3+(nelementi_foot/2*6)+6*(nnodes_foot)*(i_p-1),:)*vect_Utotf_fr_i_ep_deltaT(:,ii);
            V_Utotf_fr_i_ep_deltaT(i_p,ii)=-KKframe_fixed_foot (1+(nelementi_foot/2*6)+6*(nnodes_foot)*(i_p-1),:)*vect_Utotf_fr_i_ep_deltaT(:,ii);
        end
    end
    
    
    %DEFORMATION PARAMETERS
    for ii=1:length(Vltp_vect)
        [ X_posinflpoint_gf,DR_mag_gf,DR_mag_sag_gf(ii),DR_mag_hog_gf(ii),l_sag_gf(ii),l_hog_gf(ii) ]                    = f_reldeflection_search(1,X_foot_centr,Z_foot_centr,vect_Ugf_z                 (1+(nelementi_foot/2):(nnodes_foot):end,ii)',0.5);
        [ X_posinflpoint_fr_i,DR_mag_fr_i,DR_mag_sag_fr_i(ii),DR_mag_hog_fr_i(ii),l_sag_fr_i(ii),l_hog_fr_i(ii) ]        = f_reldeflection_search(1,X_foot_centr,Z_foot_centr,vect_Utotf_fr_i_ep_deltaT_z(1+(nelementi_foot/2):(nnodes_foot):end,ii)',0.5);
        [ X_posinflpoint_fr_i_el,DR_mag_fr_i_el,DR_mag_sag_fr_i_el(ii),DR_mag_hog_fr_i_el(ii),l_sag_fr_i_el(ii),l_hog_fr_i_el(ii) ]        = f_reldeflection_search(1,X_foot_centr,Z_foot_centr,vect_Utotf_fr_i_el_deltaT_z(1+(nelementi_foot/2):(nnodes_foot):end,ii)',0.5);
    end
    
    %Without limit
    for ii=1:length(Vltp_vect)
        [ eps_ht_gf(ii)  ,eps_hc_gf(ii)  ,pos_ht_gf(ii)  ,pos_hc_gf(ii)]            = f_hrzstr(0,X_foot_centr,Z_foot_centr,vect_Ugf_x                 (1+(nelementi_foot/2):(nnodes_foot):end,ii)); %0-1-interpolation
        [ eps_ht_fr_i(ii),eps_hc_fr_i(ii),pos_ht_fr_i(ii),pos_hc_fr_i(ii)]          = f_hrzstr(0,X_foot_centr,Z_foot_centr,vect_Utotf_fr_i_ep_deltaT_x(1+(nelementi_foot/2):(nnodes_foot):end,ii));
        [ eps_ht_fr_i_el(ii),eps_hc_fr_i_el(ii),pos_ht_fr_i_el(ii),pos_hc_fr_i_el(ii)]          = f_hrzstr(0,X_foot_centr,Z_foot_centr,vect_Utotf_fr_i_el_deltaT_x(1+(nelementi_foot/2):(nnodes_foot):end,ii));
        
    end
    
    
    %Reduction factors for the relative deflection
    % red_factor_DR_hog_free=DR_mag_hog_free./DR_mag_hog_gf;
    % red_factor_DR_sag_free=DR_mag_sag_free./DR_mag_sag_gf;
    % red_factor_DR_hog_fr_h=DR_mag_hog_fr_h./DR_mag_hog_gf;
    % red_factor_DR_sag_fr_h=DR_mag_sag_fr_h./DR_mag_sag_gf;
    red_factor_DR_hog_fr_i=DR_mag_hog_fr_i./DR_mag_hog_gf;
    red_factor_DR_sag_fr_i=DR_mag_sag_fr_i./DR_mag_sag_gf;
    
    red_factor_DR_hog_fr_i_el=DR_mag_hog_fr_i_el./DR_mag_hog_gf;
    red_factor_DR_sag_fr_i_el=DR_mag_sag_fr_i_el./DR_mag_sag_gf;
    
    %Reduction factors for the horizontal strains
    % red_factor_eht_free=eps_ht_free./eps_ht_gf;
    % red_factor_ehc_free=eps_hc_free./eps_hc_gf;
    % red_factor_eht_fr_h=eps_ht_fr_h./eps_ht_gf;
    % red_factor_ehc_fr_h=eps_hc_fr_h./eps_hc_gf;
    red_factor_eht_fr_i=eps_ht_fr_i./eps_ht_gf;
    red_factor_ehc_fr_i=eps_hc_fr_i./eps_hc_gf;
	red_factor_eht_fr_i_el=eps_ht_fr_i_el./eps_ht_gf;
    red_factor_ehc_fr_i_el=eps_hc_fr_i_el./eps_hc_gf;
    
else
    
    
    
%     %DEFORMATION PARAMETERS 

% %Without limit - using structure edges
%     for ii=1:length(Vltp_vect)
% %         [ X_posinflpoint_gf,DR_mag_gf,DR_mag_sag_gf(ii),DR_mag_hog_gf(ii),l_sag_gf(ii),l_hog_gf(ii) ]                    = f_reldeflection(0,xnod_tot,znod_tot,vect_Ugf_z                 (:,ii),0.5);
%         [ X_posinflpoint_fr_i,DR_mag_fr_i,DR_mag_sag_fr_i(ii),DR_mag_hog_fr_i(ii),l_sag_fr_i(ii),l_hog_fr_i(ii) ]        = f_reldeflection(0,xnod_tot,znod_tot,vect_Utotf_fr_i_ep_deltaT_z(:,ii),0.5);
%         [ X_posinflpoint_fr_i_el,DR_mag_fr_i_el,DR_mag_sag_fr_i_el(ii),DR_mag_hog_fr_i_el(ii),l_sag_fr_i_el(ii),l_hog_fr_i_el(ii) ]        = f_reldeflection(0,xnod_tot,znod_tot,vect_Utotf_fr_i_el_deltaT_z(:,ii),0.5);
%   
%     end
% % (Limited hogging zones) - using structure edges
%     for ii=1:length(Vltp_vect)
%         u_limit=2.5*I_gcurve_vect(ii);
%         u_limit_vect(ii)=u_limit;
%         ind_limit=(abs(xnod_tot))<u_limit;
%         [ X_posinflpoint_gf(ii,:),reldefl_mag_gf,DR_mag_sag_gf(ii),DR_mag_hog_gf(ii),l_sag_gf(ii),l_hog_gf(ii),lengthZone(ii,:) ]                    = f_reldeflection(0,xnod_tot(ind_limit),znod_tot(ind_limit),vect_Ugf_z                 ((ind_limit),ii),0.5);
% %         [ X_posinflpoint_fr_i(ii,:),DR_mag_fr_i,DR_mag_sag_fr_i(ii),DR_mag_hog_fr_i(ii),l_sag_fr_i(ii),l_hog_fr_i(ii),lengthZone(ii,:) ]        = f_reldeflection(0,xnod_tot(ind_limit),znod_tot(ind_limit),vect_Utotf_fr_i_ep_deltaT_z((ind_limit),ii),0.5);
% %         [ X_posinflpoint_fr_i_el(ii,:),DR_mag_fr_i_el,DR_mag_sag_fr_i_el(ii),DR_mag_hog_fr_i_el(ii),l_sag_fr_i_el(ii),l_hog_fr_i_el(ii),lengthZone(ii,:) ]        = f_reldeflection(0,xnod_tot(ind_limit),znod_tot(ind_limit),vect_Utotf_fr_i_el_deltaT_z((ind_limit),ii),0.5);
%     end



% %Avoid errors due to precision that results in unreal inflection points
% vect_Ugf_z(abs(vect_Ugf_z<10^-6)=0;
% Uzpg_fr_i_off_deltaT(abs(Uzpg_fr_i_off_deltaT<10^-6)=0;
% Uzpg_fr_i_off_deltaT(abs(Uzpg_fr_i_off_deltaT<10^-6)=0;
% vect_Utotf_fr_i_ep_deltaT_z(abs(vect_Utotf_fr_i_ep_deltaT_z<10^-6)=0;
% vect_Utotf_fr_i_el_deltaT_z(abs(vect_Utotf_fr_i_el_deltaT_z<10^-6)=0;



    for ii=1:length(Vltp_vect)
        [ X_posinflpoint_gf     ,DR_mag_gf     ,DR_mag_sag_gf(ii)     ,DR_mag_hog_gf(ii)     ,l_sag_gf(ii)     ,l_hog_gf(ii)     ,lengthZone_gf     ,X_posinflpoint_nosearch_gf     ,X_posinflpoint_searched_gf]          	= f_reldeflection_search(0,xnod_tot,znod_tot,vect_Ugf_z                 (:,ii),0.5);
        
        if switch_NA==1
        [ X_posinflpoint_fr_i   ,DR_mag_fr_i   ,DR_mag_sag_fr_i(ii)   ,DR_mag_hog_fr_i(ii)   ,l_sag_fr_i(ii)   ,l_hog_fr_i(ii)   ,lengthZone_fr_i   ,X_posinflpoint_nosearch_fr_i   ,X_posinflpoint_searched_fr_i   ]    	= f_reldeflection_search(0,xnod_tot,znod_tot,Uzpg_fr_i_off_deltaT(:,ii),0.5);
        [ X_posinflpoint_fr_i_el,DR_mag_fr_i_el,DR_mag_sag_fr_i_el(ii),DR_mag_hog_fr_i_el(ii),l_sag_fr_i_el(ii),l_hog_fr_i_el(ii),lengthZone_fr_i_el,X_posinflpoint_nosearch_fr_i_el,X_posinflpoint_searched_fr_i_el ]   	= f_reldeflection_search(0,xnod_tot,znod_tot,Uzpg_fr_i_off_deltaT(:,ii),0.5);
        else
        [ X_posinflpoint_fr_i   ,DR_mag_fr_i   ,DR_mag_sag_fr_i(ii)   ,DR_mag_hog_fr_i(ii)   ,l_sag_fr_i(ii)   ,l_hog_fr_i(ii)   ,lengthZone_fr_i   ,X_posinflpoint_nosearch_fr_i   ,X_posinflpoint_searched_fr_i   ]    	= f_reldeflection_search(0,xnod_tot,znod_tot,vect_Utotf_fr_i_ep_deltaT_z(:,ii),0.5);
        [ X_posinflpoint_fr_i_el,DR_mag_fr_i_el,DR_mag_sag_fr_i_el(ii),DR_mag_hog_fr_i_el(ii),l_sag_fr_i_el(ii),l_hog_fr_i_el(ii),lengthZone_fr_i_el,X_posinflpoint_nosearch_fr_i_el,X_posinflpoint_searched_fr_i_el ]   	= f_reldeflection_search(0,xnod_tot,znod_tot,vect_Utotf_fr_i_el_deltaT_z(:,ii),0.5);
        end

        if l_sag_fr_i(ii)>0
            ind_sag=find(DR_mag_fr_i==DR_mag_sag_fr_i(ii));
            ind_sag=ind_sag(1);
            if     X_posinflpoint_fr_i([ind_sag])./I_gcurve_vect(ii)<-1 && X_posinflpoint_fr_i([ind_sag+1])./I_gcurve_vect(ii)<-1
                ind_correction_sag(ii)=0;
            elseif X_posinflpoint_fr_i([ind_sag])./I_gcurve_vect(ii)>+1 && X_posinflpoint_fr_i([ind_sag+1])./I_gcurve_vect(ii)>+1
                ind_correction_sag(ii)=0;
            else
                ind_correction_sag(ii)=1;
            end
        else
            ind_correction_sag(ii)=1; %to not consider sagging far away from tunnel
        end
                
        if l_sag_fr_i_el(ii)>0
            ind_sag_el=find(DR_mag_fr_i_el==DR_mag_sag_fr_i_el(ii));
            ind_sag_el=ind_sag_el(1);
            if  X_posinflpoint_fr_i_el(ind_sag_el)./I_gcurve_vect(ii)<-1 && X_posinflpoint_fr_i_el(ind_sag_el+1)./I_gcurve_vect(ii)<-1
                ind_correction_sag_el(ii)=0;
            elseif X_posinflpoint_fr_i_el(ind_sag_el)./I_gcurve_vect(ii)>+1 && X_posinflpoint_fr_i_el(ind_sag_el+1)./I_gcurve_vect(ii)>+1
                ind_correction_sag_el(ii)=0;
            else
                ind_correction_sag_el(ii)=1;
            end
        else
            ind_correction_sag_el(ii)=1;
        end
        
        
        if num_storey>0
        else
          
            
        if switch_NA==1
        [ eps_avg_fr_i_gf,eps_avg_sag_gf(ii)     ,eps_avg_hog_gf(ii)    ]  =    f_hrzstr_avg_DR(xnod_tot',znod_tot', vect_Ugf_x                 (:,ii)',X_posinflpoint_searched_gf     ,DR_mag_gf);
        [ eps_avg_fr_i   ,eps_avg_sag_fr_i(ii)   ,eps_avg_hog_fr_i(ii)   ] =    f_hrzstr_avg_DR(xnod_tot',znod_tot',Uxpg_fr_i_off_deltaT(:,ii)',X_posinflpoint_searched_fr_i   ,DR_mag_fr_i);  
        [ eps_avg_fr_i_el,eps_avg_sag_fr_i_el(ii),eps_avg_hog_fr_i_el(ii)] =    f_hrzstr_avg_DR(xnod_tot',znod_tot',Uxpg_fr_i_off_deltaT(:,ii)',X_posinflpoint_searched_fr_i_el,DR_mag_fr_i_el);
        else
        [ eps_avg_fr_i_gf,eps_avg_sag_gf(ii)     ,eps_avg_hog_gf(ii)    ]  =    f_hrzstr_avg_DR(xnod_tot',znod_tot', vect_Ugf_x                 (:,ii)',X_posinflpoint_searched_gf     ,DR_mag_gf);
        [ eps_avg_fr_i   ,eps_avg_sag_fr_i(ii)   ,eps_avg_hog_fr_i(ii)   ] =    f_hrzstr_avg_DR(xnod_tot',znod_tot',vect_Utotf_fr_i_ep_deltaT_x(:,ii)',X_posinflpoint_searched_fr_i   ,DR_mag_fr_i);  
        [ eps_avg_fr_i_el,eps_avg_sag_fr_i_el(ii),eps_avg_hog_fr_i_el(ii)] =    f_hrzstr_avg_DR(xnod_tot',znod_tot',vect_Utotf_fr_i_el_deltaT_x(:,ii)',X_posinflpoint_searched_fr_i_el,DR_mag_fr_i_el);
        end
            
        end
    
    end
    
    
    for ii=1:length(Vltp_vect)
        
        if switch_NA==1
        [ eps_ht_gf(ii)  ,eps_hc_gf(ii)  ,pos_ht_gf(ii)  ,pos_hc_gf(ii)]                        = f_hrzstr(0,xnod_tot,znod_tot,vect_Ugf_x                 (:,ii)); %0-1-interpolation
        [ eps_ht_fr_i(ii),eps_hc_fr_i(ii),pos_ht_fr_i(ii),pos_hc_fr_i(ii)]                      = f_hrzstr(0,xnod_tot,znod_tot,Uxpg_fr_i_off_deltaT(:,ii));
        [ eps_ht_fr_i_el(ii),eps_hc_fr_i_el(ii),pos_ht_fr_i_el(ii),pos_hc_fr_i_el(ii)]          = f_hrzstr(0,xnod_tot,znod_tot,Uxpg_fr_i_off_deltaT(:,ii));
        else
        [ eps_ht_gf(ii)  ,eps_hc_gf(ii)  ,pos_ht_gf(ii)  ,pos_hc_gf(ii)]                        = f_hrzstr(0,xnod_tot,znod_tot,vect_Ugf_x                 (:,ii)); %0-1-interpolation
        [ eps_ht_fr_i(ii),eps_hc_fr_i(ii),pos_ht_fr_i(ii),pos_hc_fr_i(ii)]                      = f_hrzstr(0,xnod_tot,znod_tot,vect_Utotf_fr_i_ep_deltaT_x(:,ii));
        [ eps_ht_fr_i_el(ii),eps_hc_fr_i_el(ii),pos_ht_fr_i_el(ii),pos_hc_fr_i_el(ii)]          = f_hrzstr(0,xnod_tot,znod_tot,vect_Utotf_fr_i_el_deltaT_x(:,ii));
        end
        
        
 end
    
    
 
%     %Reduction factors for the relative deflection
%     red_factor_DR_hog_fr_i=DR_mag_hog_fr_i./DR_mag_hog_gf;
%     red_factor_DR_sag_fr_i=DR_mag_sag_fr_i./DR_mag_sag_gf;
%     red_factor_DR_hog_fr_i_el=DR_mag_hog_fr_i_el./DR_mag_hog_gf;
%     red_factor_DR_sag_fr_i_el=DR_mag_sag_fr_i_el./DR_mag_sag_gf; 
%     %Reduction factors for the horizontal strains
%     red_factor_eht_fr_i=eps_ht_fr_i./eps_ht_gf;
%     red_factor_ehc_fr_i=eps_hc_fr_i./eps_hc_gf;
%     red_factor_eht_fr_i_el=eps_ht_fr_i_el./eps_ht_gf;
%     red_factor_ehc_fr_i_el=eps_hc_fr_i_el./eps_hc_gf;  

    %Reduction factors for the relative deflection
    red_factor_DR_hog_fr_i=DR_mag_hog_fr_i./DR_mag_hog_gf;
    red_factor_DR_sag_fr_i=ind_correction_sag.*  DR_mag_sag_fr_i./DR_mag_sag_gf;
    red_factor_DR_hog_fr_i_el=DR_mag_hog_fr_i_el./DR_mag_hog_gf;
    red_factor_DR_sag_fr_i_el=ind_correction_sag_el.* DR_mag_sag_fr_i_el./DR_mag_sag_gf; 
    l_sag_fr_i   =ind_correction_sag   .*l_sag_fr_i   ;
    l_sag_fr_i_el=ind_correction_sag_el.*l_sag_fr_i_el;
    %Reduction factors for the horizontal strains
    red_factor_eht_fr_i=eps_ht_fr_i./eps_ht_gf;
    red_factor_ehc_fr_i= eps_hc_fr_i./eps_hc_gf;
    red_factor_eht_fr_i_el=eps_ht_fr_i_el./eps_ht_gf;
    red_factor_ehc_fr_i_el= eps_hc_fr_i_el./eps_hc_gf;  

        if num_storey>0
        else
        red_factor_eh_avg_sag_fr_i    =ind_correction_sag    .*eps_avg_sag_fr_i   ./eps_avg_sag_gf;
        red_factor_eh_avg_hog_fr_i    =eps_avg_hog_fr_i   ./eps_avg_hog_gf;
        red_factor_eh_avg_sag_fr_i_el =ind_correction_sag_el .*eps_avg_sag_fr_i_el./eps_avg_sag_gf;
        red_factor_eh_avg_hog_fr_i_el =eps_avg_hog_fr_i_el./eps_avg_hog_gf;  
        end
    
    DRsag_gauss=0.24657034*Vltp_vect/100*(2*Rt/ht)^2;
	N_DR_hog_gf=DR_mag_hog_gf./DRsag_gauss;
    N_DR_sag_gf=DR_mag_sag_gf./DRsag_gauss;
	N_DR_hog_fr_i=DR_mag_hog_fr_i./DRsag_gauss;
    N_DR_sag_fr_i=DR_mag_sag_fr_i./DRsag_gauss;
    N_DR_hog_fr_i_el=DR_mag_hog_fr_i_el./DRsag_gauss;
    N_DR_sag_fr_i_el=DR_mag_sag_fr_i_el./DRsag_gauss; 
    
    ro_sag= Efoot.*(dfoot.^3./12)./Es./l_sag_gf.^3;   
    ro_hog= Efoot.*(dfoot.^3./12)./Es./l_hog_gf.^3;   

    B_i_lt=min(xnod_tot)./I_gcurve_vect;
    B_i_rt=max(xnod_tot)./I_gcurve_vect;
    
    
    Br_temp=B_i_rt;
    Bl_temp=B_i_lt;
    B_i_max=max([abs(Br_temp),abs(Bl_temp)],[],2);
    B_i_min= sign(Br_temp./Bl_temp) .*min([abs(Br_temp),abs(Bl_temp)],[],2);
    

    EI_foot=Efoot*bfoot*dfoot^3/12;
    EA_foot=Efoot*bfoot*dfoot;
    
    e_B=X_foot_centr/B_X_foot;
    
    k_foot=10*(1+ni_foot)/(12+11*ni_foot);
    
    if max(max(F_soil_fr_i_el_T_z))>0
    ind_gap_EL=1;
    else
    ind_gap_EL=0;
    end
    
    
    result_TSI=[... %result_red_factor
        EI_foot/Es/I_gcurve_vect(end)^3,...
        EI_foot,...
        EA_foot,...
        nis,...
        e_B,... %5
        B_i_lt(end),...
        B_i_rt(end),...
        B_i_max(end),...
        B_i_min(end),...
        qz_foot/Es/bfoot,... %10
        l_sag_gf(end)/I_gcurve_vect(end),...
        l_hog_gf(end)/I_gcurve_vect(end),...
        l_sag_fr_i_el(end)/I_gcurve_vect(end),... 
        l_hog_fr_i_el(end)/I_gcurve_vect(end),...
        ro_sag(end),...%15
        ro_hog(end), ...
        red_factor_DR_sag_fr_i_el(end),...
        red_factor_DR_hog_fr_i_el(end),... 
        N_DR_sag_fr_i_el(end),...
        N_DR_hog_fr_i_el(end),...%20
        I_gcurve_vect(end),...
        l_sag_fr_i_el(end)/l_sag_gf(end),... % 22 MBsab
        l_hog_fr_i_el(end)/l_hog_gf(end),... %23 MBhog     
        EGratio/(2*(1+ni_foot)),...
        coeff_TIM*12*EI_foot/(B_X_foot(1)^2*k_foot*EA_foot/EGratio),...  % 25 F=12kEI/2(1+ni)GAsh^2
        B_X_foot(1)./bfoot,...
        coeff_TIM*12*EI_foot/(l_sag_gf(end)^2*k_foot*EA_foot/EGratio),...
        coeff_TIM*12*EI_foot/(l_hog_gf(end)^2*k_foot*EA_foot/EGratio),...
        ind_gap_EL,...
        max(max(F_soil_fr_i_el_T_z)),... % 30 
        Rt/B_X_foot(1),...
        red_factor_DR_sag_fr_i(end),...
        red_factor_DR_hog_fr_i(end),...
        EI_foot/bfoot    /Es/B_X_foot(1)^3,... %Values for the skew beam % 34
        EA_foot          /Es/B_X_foot(1)^1,...
        k_foot,...
        ni_foot,...
        c_shear,...
        s_shear,...
        B_X_foot(1)./bfoot,...% 40
        H_target/B_X_foot(1),...
        d_na/H_target,...
        EGratio,...
        ht/B_X_foot(1),...
        theta_degree,...% 45
        dorig/B_X_foot(1),...
        Kx,...
        Ky,...
        delta_ratio,...
        ];  


    %SAVE
    if exist('Data_TSI.mat')==0
        data_TSI(1,:)=result_TSI;
        save('Data_TSI.mat','data_TSI')
    else
        load Data_TSI.mat
        i=1+size(data_TSI,1);
        data_TSI(i,:)=result_TSI;
        save('Data_TSI.mat','data_TSI')
    end
    
end









%%





%% Calculation of internal forces
% 

if nfoot==1
% % % for ii=1:nnodes_foot-1
% % % 
% % % Fin_P=KBern3Delt*vect_Utotf_fr_i_el_P(1+6*(ii-1):12+6*(ii-1),end);
% % % 
% % % Fin_P_1(1:6,ii)=Fin_P(1:6);
% % % Fin_P_2(1:6,ii)=Fin_P(7:12);
% % % 
% % % Fin_T=KBern3Delt*vect_Utotf_fr_i_el_T(1+6*(ii-1):12+6*(ii-1),end);
% % % 
% % % Fin_T_1(1:6,ii)=Fin_T(1:6);
% % % Fin_T_2(1:6,ii)=Fin_T(7:12);
% % % 
% % % Fin_deltaT_1=Fin_T_1-Fin_P_1;
% % % Fin_deltaT_2=Fin_T_2-Fin_P_2;
% % % 
% % % end
% % % 
% % % Xord(1:2:2*(nnodes_foot-1))=xnod_tot(1:end-1);
% % % Xord(2:2:2*(nnodes_foot-1))=xnod_tot(2:end);
% % % 
% % % %For the footing/equivalent beam I have adopted the traditional reference
% % % %system for 
% % % 
% % % F_M_deltaT_el(1:2:2*(nnodes_foot-1))=-Fin_deltaT_1(5,:);
% % % F_M_deltaT_el(2:2:2*(nnodes_foot-1))=Fin_deltaT_2(5,:);
% % % F_N_deltaT_el(1:2:2*(nnodes_foot-1))=-Fin_deltaT_1(1,:);
% % % F_N_deltaT_el(2:2:2*(nnodes_foot-1))=Fin_deltaT_2(1,:);
% % % F_S_deltaT_el(1:2:2*(nnodes_foot-1))=-Fin_deltaT_1(3,:);
% % % F_S_deltaT_el(2:2:2*(nnodes_foot-1))=Fin_deltaT_2(3,:);
% % % 
% % % 
% % % F_M_P_el(1:2:2*(nnodes_foot-1))=-Fin_P_1(5,:);
% % % F_M_P_el(2:2:2*(nnodes_foot-1))=Fin_P_2(5,:);
% % % F_N_P_el(1:2:2*(nnodes_foot-1))=-Fin_P_1(1,:);
% % % F_N_P_el(2:2:2*(nnodes_foot-1))=Fin_P_2(1,:);
% % % F_S_P_el(1:2:2*(nnodes_foot-1))=-Fin_P_1(3,:);
% % % F_S_P_el(2:2:2*(nnodes_foot-1))=Fin_P_2(3,:);
for jj=1:num_incr_P
for ii=1:nnodes_foot-1

    
if switch_NA==1
    Fin_P=KBern3Delt*U_node_fr_i_off_P(1+6*(ii-1):12+6*(ii-1),jj);
else
    Fin_P=KBern3Delt*vect_Utotf_fr_i_el_P(1+6*(ii-1):12+6*(ii-1),jj);
end


Fin_P_1(1:6,ii)=Fin_P(1:6);
Fin_P_2(1:6,ii)=Fin_P(7:12);


end

Xord(1:2:2*(nnodes_foot-1))=xnod_tot(1:end-1);
Xord(2:2:2*(nnodes_foot-1))=xnod_tot(2:end);

%For the footing/equivalent beam I have adopted the traditional reference
%system for 

F_M_P_el(1:2:2*(nnodes_foot-1))=-Fin_P_1(5,:);
F_M_P_el(2:2:2*(nnodes_foot-1))=Fin_P_2(5,:);
F_N_P_el(1:2:2*(nnodes_foot-1))=-Fin_P_1(1,:);
F_N_P_el(2:2:2*(nnodes_foot-1))=Fin_P_2(1,:);
F_S_P_el(1:2:2*(nnodes_foot-1))=-Fin_P_1(3,:);
F_S_P_el(2:2:2*(nnodes_foot-1))=Fin_P_2(3,:);

F_M_P_el_M(:,jj)=F_M_P_el;
F_N_P_el_M(:,jj)=F_N_P_el;
F_S_P_el_M(:,jj)=F_S_P_el;
end

for jj=1:num_incr_TUNNEL
for ii=1:nnodes_foot-1

if switch_NA==1
    Fin_P=KBern3Delt*U_node_fr_i_off_P(1+6*(ii-1):12+6*(ii-1),end);
else
    Fin_P=KBern3Delt*vect_Utotf_fr_i_el_P(1+6*(ii-1):12+6*(ii-1),end);
end
    


Fin_P_1(1:6,ii)=Fin_P(1:6);
Fin_P_2(1:6,ii)=Fin_P(7:12);

if switch_NA==1
    Fin_T=KBern3Delt*U_node_fr_i_off(1+6*(ii-1):12+6*(ii-1),jj);
else
    Fin_T=KBern3Delt*vect_Utotf_fr_i_el_T(1+6*(ii-1):12+6*(ii-1),jj);
end

Fin_T_1(1:6,ii)=Fin_T(1:6);
Fin_T_2(1:6,ii)=Fin_T(7:12);



end

Fin_deltaT_1=Fin_T_1-Fin_P_1;
Fin_deltaT_2=Fin_T_2-Fin_P_2;

Xord(1:2:2*(nnodes_foot-1))=xnod_tot(1:end-1);
Xord(2:2:2*(nnodes_foot-1))=xnod_tot(2:end);

%For the footing/equivalent beam I have adopted the traditional reference
%system for 

F_M_deltaT_el(1:2:2*(nnodes_foot-1))=-Fin_deltaT_1(5,:);
F_M_deltaT_el(2:2:2*(nnodes_foot-1))=Fin_deltaT_2(5,:);
F_N_deltaT_el(1:2:2*(nnodes_foot-1))=-Fin_deltaT_1(1,:);
F_N_deltaT_el(2:2:2*(nnodes_foot-1))=Fin_deltaT_2(1,:);
F_S_deltaT_el(1:2:2*(nnodes_foot-1))=-Fin_deltaT_1(3,:);
F_S_deltaT_el(2:2:2*(nnodes_foot-1))=Fin_deltaT_2(3,:);

F_M_deltaT_el_M(:,jj)=F_M_deltaT_el;
F_N_deltaT_el_M(:,jj)=F_N_deltaT_el;
F_S_deltaT_el_M(:,jj)=F_S_deltaT_el;

end

if switch_plot==0
else
    
figure(22);hold on;
h(1)=plot(Xord,F_M_deltaT_el,'-k');
h(2)=plot(Xord,F_N_deltaT_el,'-r');
h(3)=plot(Xord,F_S_deltaT_el,'-b');

h(1)=plot(Xord,F_M_deltaT_el_M(:,end),'-ok');
h(2)=plot(Xord,F_N_deltaT_el_M(:,end),'-or');
h(3)=plot(Xord,F_S_deltaT_el_M(:,end),'-ob');

set(gca,'Ydir','reverse')


Legend1 =legend( [h(1),h(2),h(3)],...
    {'M','N','V'},...
    'Location','southeastoutside');%,'PlotBoxAspectRatio',ratio_leg21'Position',pos_leg_1','PlotBoxAspectRatioMode','manual'
legend('boxon')
ylabel('Forces (N,Nm)');% 
xlabel('Offset (m)');%

title('Internal forces for the elastic solution')
end
%% Bending and shear deformations


if switch_NA==1
    utot=U_node_fr_i_off_deltaT;
else
    utot=vect_Utotf_fr_i_el_deltaT;
end

% % % % 
% % % % I_fot=(bfoot.*dfoot.^3./12);
% % % % du_bend=(-(F_M_deltaT_el_M(1:2:end,:)/2+F_M_deltaT_el_M(2:2:end,:)/2)./(Efoot.*I_fot)).*(ones(nnodes_foot-1,num_incr_TUNNEL)*h_el_foot^2/2);
% % % % du_shear=(F_S_deltaT_el_M(1:2:end,:)./((10+ni_foot*10)./(12+11*ni_foot)*bfoot*dfoot)./(Efoot./EGratio)).*(ones(nnodes_foot-1,num_incr_TUNNEL)*h_el_foot);
% % % % du_tilt=-utot(5:6:end-6,:).*(ones(nnodes_foot-1,num_incr_TUNNEL)*h_el_foot);
% % % % 
% % % % u_bend  =cumsum(du_bend(:,:),1);
% % % % u_shear =cumsum(du_shear(:,:),1);%
% % % % u_tilt  =cumsum(du_tilt(:,:),1);%
% % % % 
% % % % u_bend =[zeros(1,num_incr_TUNNEL);u_bend];
% % % % u_shear=[zeros(1,num_incr_TUNNEL);u_shear];
% % % % u_tilt =[zeros(1,num_incr_TUNNEL);u_tilt];
% % % % 
% % % % curvature=du_bend./(h_el_foot^2/2);
% % % % gamma_strain=du_shear./(h_el_foot);



%Deformations from forces
I_fot=(bfoot.*dfoot.^3./12);
curvature=(-(F_M_deltaT_el_M(1:2:end,:)/2+F_M_deltaT_el_M(2:2:end,:)/2)./(Efoot.*I_fot));
gamma_strain=(F_S_deltaT_el_M(1:2:end,:)./((10+ni_foot*10)./(12+11*ni_foot)*bfoot*dfoot)./(Efoot./EGratio));

du_bend=curvature.*(ones(nnodes_foot-1,num_incr_TUNNEL)*h_el_foot^2/2);
du_shear=gamma_strain.*(ones(nnodes_foot-1,num_incr_TUNNEL)*h_el_foot);
du_tilt=-utot(5:6:end-6,:).*(ones(nnodes_foot-1,num_incr_TUNNEL)*h_el_foot);

u_bend  =cumsum(du_bend(:,:),1);
u_shear =cumsum(du_shear(:,:),1);%
u_tilt  =cumsum(du_tilt(:,:),1);%

u_bend =[zeros(1,num_incr_TUNNEL);u_bend];
u_shear=[zeros(1,num_incr_TUNNEL);u_shear];
u_tilt =[zeros(1,num_incr_TUNNEL);u_tilt];





%Tilt connecting edges
lbuild=(max(x_global_foot(end))-min(x_global_foot(1)));
tetha_tilt=atan((utot(end-3,:)-utot(3,:))/lbuild);
du_tiltbody = repmat(tetha_tilt,nnodes_foot-1,1)*h_el_foot;%
u_tiltbody  =cumsum(du_tiltbody(:,:),1);%
u_tiltbody =[zeros(1,num_incr_TUNNEL);u_tiltbody];


%Structure fully flexible at shear
du_tot=diff(utot(3:6:end,:),1);
du_shear_fullflexshear=du_tot-du_tiltbody;
gamma_strain_fullflexshear=du_shear_fullflexshear./(h_el_foot);

%Deformations from displacement vector
du_bend2=-(utot(5+6:6:end,:)-utot(5:6:end-6,:)).*(ones(nnodes_foot-1,num_incr_TUNNEL)*h_el_foot/2);
u_bend2  =cumsum(du_bend2(:,:),1);
u_bend2 =[zeros(1,num_incr_TUNNEL);u_bend2];
u_shear2=utot(3:6:end,:)-repmat(utot(3,:),nnodes_foot,1)-u_bend2-u_tilt;
du_shear2=diff(u_shear2,1,1);

omegaaa=utot(5:6:end,:);

u_bend_bernoulli=utot(3:6:end,:)-repmat(utot(3,:),nnodes_foot,1)-u_tilt;
du_bend_bernoulli=diff(u_bend_bernoulli,1,1);
curvature_bernoulli=du_bend_bernoulli./(h_el_foot^2/2);



if switch_plot==0
else
    
    
figure(219);hold on;
title('Total Shear and Bending Deflections')
h(1)=plot(xnod_tot,u_bend(:,end),'-k');
h(2)=plot(xnod_tot,u_shear(:,end),'-r');
h(3)=plot(xnod_tot,u_tilt(:,end),'-b');
% h(3)=plot(xnod_tot,u_tiltbody(:,end),'--b');
% h(3)=plot(xnod_tot,u_tilt_bending(:,end),'-xb');
% h(3)=plot(xnod_tot,u_tilt_first(:,end),'-xr');
% h(3)=plot(xnod_tot,u_tilt_bendingandfirst(:,end),'-xk');

h(4)=plot(xnod_tot,u_bend(:,end)+u_shear(:,end)+u_tilt(:,end),'-c');
h(5)=plot(xnod_tot,utot(3:6:end,end)-utot(3,end),'og');
h(5)=plot(xnod_tot,utot(3:6:end,end),'og');
h(5)=plot(xnod_tot,0*utot(3:6:end,end),'k:');

set(gca,'Ydir','reverse')





figure(250);hold on;
h(1)=plot(xnod_tot,u_bend(:,end),'ok');
h(2)=plot(xnod_tot,u_bend2(:,end),'-k');
h(1)=plot(xnod_tot,u_shear(:,end),'or');
h(2)=plot(xnod_tot,u_shear2(:,end),'-r');


figure(251);hold on;
h(1)=plot(xnod_tot(1:1:end-1)/2+xnod_tot(2:1:end)/2,du_bend(:,end),'ok');
h(2)=plot(xnod_tot(1:1:end-1)/2+xnod_tot(2:1:end)/2,du_bend2(:,end),'-k');
h(1)=plot(xnod_tot(1:1:end-1)/2+xnod_tot(2:1:end)/2,du_shear(:,end),'or');
h(2)=plot(xnod_tot(1:1:end-1)/2+xnod_tot(2:1:end)/2,du_shear2(:,end),'-r');

h(3)=plot(xnod_tot(1:1:end-1)/2+xnod_tot(2:1:end)/2,du_shear_fullflexshear(:,end),'-c');


figure(301);
subplot(1,2,1);hold on;
plot(xnod_tot(1:1:end-1)/2+xnod_tot(2:1:end)/2,curvature(:,end),'k');
% plot(xnod_tot(1:1:end-1)/2+xnod_tot(2:1:end)/2,curvature_bernoulli(:,end),':k');
subplot(1,2,2)
plot(xnod_tot(1:1:end-1)/2+xnod_tot(2:1:end)/2,gamma_strain(:,end),'r')
set(gca,'Ydir','reverse')

figure(239);hold on;
title('Element Deflections')
h(1)=plot(xnod_tot(1:1:end-1)/2+xnod_tot(2:1:end)/2,du_bend(:,end),'ok');
% h(10)=plot(xnod_tot(1:1:end-1)/2+xnod_tot(2:1:end)/2,du_bend_bernoulli(:,end),'-ok');
h(2)=plot(xnod_tot(1:1:end-1)/2+xnod_tot(2:1:end)/2,du_shear(:,end),'or');
h(3)=plot(xnod_tot(1:1:end-1)/2+xnod_tot(2:1:end)/2,du_tilt(:,end),'ob');


Legend1 =legend( [h(1),h(2),h(3)],...
    {'Bending','Shearing','Tilting'},...
    'Location','southeastoutside');%,'PlotBoxAspectRatio',ratio_leg21'Position',pos_leg_1','PlotBoxAspectRatioMode','manual'
legend('boxon')
end


%% Strain assessment for the equivalent beam for the ELASTIC SOLUTION

%($t$ and $s$ have positive signs below the neutral axis while tensile fibres 
%are at the beam bottom for positive $M$)

I_fot=(bfoot.*dfoot.^3./12);
% t_fot=(dfoot/2);

% epsilon_dmax=(c_shear).*abs(F_S_deltaT_el)./((10+ni_foot*10)./(12+11*ni_foot)*bfoot*dfoot)./(2*Efoot./EGratio);
% epsilon_bmax_top=F_M_deltaT_el.*(-1*(H_target-d_na))/(Efoot.*I_fot);
% epsilon_bmax_bot=F_M_deltaT_el.*(d_na)/(Efoot.*I_fot);
% epsilon_bmax_pri=abs(F_M_deltaT_el).*(s_shear)/(Efoot.*I_fot);
% % epsilon_bmax=abs(F_M_deltaT_el).*(dfoot/2)/(Efoot.*I_fot);
% epsilon_h=F_N_deltaT_el./(Efoot.*(bfoot.*dfoot));
 

epsilon_dmax=(c_shear).*abs(F_S_deltaT_el_M(:,:))./((10+ni_foot*10)./(12+11*ni_foot)*bfoot*dfoot)./(2*Efoot./EGratio);
epsilon_bmax_top=F_M_deltaT_el_M(:,:).*(-1*(H_target-d_na))/(Efoot.*I_fot);
epsilon_bmax_bot=F_M_deltaT_el_M(:,:).*(d_na)/(Efoot.*I_fot);
epsilon_bmax_pri=abs(F_M_deltaT_el_M(:,:)).*(s_shear)/(Efoot.*I_fot);
% epsilon_bmax=abs(F_M_deltaT_el).*(dfoot/2)/(Efoot.*I_fot);
epsilon_h=F_N_deltaT_el_M(:,:)./(Efoot.*(bfoot.*dfoot));

% epsilon_br=epsilon_bmax+epsilon_h;
% epsilon_dr=epsilon_h.*((1-ni_foot)/2)+(epsilon_h.^2*((1+ni_foot)/2)^2+epsilon_dmax.^2).^0.5;   

epsilon_br_top=epsilon_bmax_top+epsilon_h;
epsilon_br_bot=epsilon_bmax_bot+epsilon_h;
epsilon_br_env=max(epsilon_br_top,epsilon_br_bot);

epsilon_dr=(epsilon_bmax_pri+epsilon_h).*((1-ni_foot)/2)+((epsilon_bmax_pri+epsilon_h).^2*((1+ni_foot)/2)^2+epsilon_dmax.^2).^0.5;   

epsilon_brmax_fba=max(epsilon_br_env);
epsilon_drmax_fba=max(epsilon_dr);

epsion_max_fba=max([epsilon_brmax_fba;epsilon_drmax_fba]);


%DISPLACEMENT BASED APPROACH
I_fot_1=(1.*dfoot.^3./12); 
t_fot=(d_na);
epsilon_bmax_bur_sag_gf=-DR_mag_sag_gf     (:)./(l_sag_gf(:)     ./12./t_fot+3.*I_fot_1.*EGratio./2./t_fot./l_sag_gf(:)     ./H_target);
epsilon_bmax_bur_sag_bl=-DR_mag_sag_fr_i_el(:)./(l_sag_fr_i_el(:)./12./t_fot+3.*I_fot_1.*EGratio./2./t_fot./l_sag_fr_i_el(:)./H_target);
epsilon_dmax_bur_sag_gf=-DR_mag_sag_gf     (:)./(1+H_target.*l_sag_gf(:)     .^2./18./I_fot_1./EGratio);
epsilon_dmax_bur_sag_bl=-DR_mag_sag_fr_i_el(:)./(1+H_target.*l_sag_fr_i_el(:).^2./18./I_fot_1./EGratio);

I_fot_1=(1.*dfoot.^3./12);
t_fot=(H_target-d_na);
epsilon_bmax_bur_hog_gf=DR_mag_hog_gf      (:)./(l_hog_gf(:)     ./12./t_fot+3.*I_fot_1.*EGratio./2./t_fot./l_hog_gf(:)     ./H_target);
epsilon_bmax_bur_hog_bl=DR_mag_hog_fr_i_el (:)./(l_hog_fr_i_el(:)./12./t_fot+3.*I_fot_1.*EGratio./2./t_fot./l_hog_fr_i_el(:)./H_target);
epsilon_dmax_bur_hog_gf=DR_mag_hog_gf      (:)./(1+H_target.*l_hog_gf(:)     .^2./18./I_fot_1./EGratio);
epsilon_dmax_bur_hog_bl=DR_mag_hog_fr_i_el (:)./(1+H_target.*l_hog_fr_i_el(:).^2./18./I_fot_1./EGratio);
 
if num_storey==0
Epsilon_h=eps_avg_hog_fr_i_el(:);
Epsilon_bmax=epsilon_bmax_bur_hog_bl;
Epsilon_dmax=epsilon_dmax_bur_hog_bl;
epsilon_br_dba_hog_bl=Epsilon_bmax+Epsilon_h;
epsilon_dr_dba_hog_bl=Epsilon_h.*((1-ni_foot)./2)+(Epsilon_h.^2.*((1+ni_foot)./2).^2+Epsilon_dmax.^2).^0.5;   

Epsilon_h=eps_avg_sag_fr_i_el(:);
Epsilon_bmax=epsilon_bmax_bur_sag_bl;
Epsilon_dmax=epsilon_dmax_bur_sag_bl;
epsilon_br_dba_sag_bl=Epsilon_bmax+Epsilon_h;
epsilon_dr_dba_sag_bl=Epsilon_h.*((1-ni_foot)./2)+(Epsilon_h.^2.*((1+ni_foot)./2).^2+Epsilon_dmax.^2).^0.5; 

Epsilon_h=eps_avg_hog_gf(:);
Epsilon_bmax=epsilon_bmax_bur_hog_gf;
Epsilon_dmax=epsilon_dmax_bur_hog_gf;
epsilon_br_dba_hog_gf=Epsilon_bmax+Epsilon_h;
epsilon_dr_dba_hog_gf=Epsilon_h.*((1-ni_foot)./2)+(Epsilon_h.^2.*((1+ni_foot)./2).^2+Epsilon_dmax.^2).^0.5;   

Epsilon_h=eps_avg_sag_gf(:);
Epsilon_bmax=epsilon_bmax_bur_sag_gf;
Epsilon_dmax=epsilon_dmax_bur_sag_gf;
epsilon_br_dba_sag_gf=Epsilon_bmax+Epsilon_h;
epsilon_dr_dba_sag_gf=Epsilon_h.*((1-ni_foot)./2)+(Epsilon_h.^2.*((1+ni_foot)./2).^2+Epsilon_dmax.^2).^0.5; 




%MODIFICATION FACTOR BASED APPROACH
I_fot_1=(1.*dfoot.^3./12); 
t_fot=(d_na);
epsilon_bmax_mba_sag_gf=-DR_mag_sag_gf     (:)./(l_sag_gf(:)     ./12./t_fot+3.*I_fot_1.*EGratio./2./t_fot./l_sag_gf(:)     ./H_target);
epsilon_bmax_mba_sag_bl=-DR_mag_sag_fr_i_el(:)./(l_sag_gf(:)     ./12./t_fot+3.*I_fot_1.*EGratio./2./t_fot./l_sag_gf(:)./H_target);
epsilon_dmax_mba_sag_gf=-DR_mag_sag_gf     (:)./(1+H_target.*l_sag_gf(:)     .^2./18./I_fot_1./EGratio);
epsilon_dmax_mba_sag_bl=-DR_mag_sag_fr_i_el(:)./(1+H_target.*l_sag_gf(:)     .^2./18./I_fot_1./EGratio);

I_fot_1=(1.*dfoot.^3./12);
t_fot=(H_target-d_na);
epsilon_bmax_mba_hog_gf=DR_mag_hog_gf      (:)./(l_hog_gf(:)     ./12./t_fot+3.*I_fot_1.*EGratio./2./t_fot./l_hog_gf(:)     ./H_target);
epsilon_bmax_mba_hog_bl=DR_mag_hog_fr_i_el (:)./(l_hog_gf(:)     ./12./t_fot+3.*I_fot_1.*EGratio./2./t_fot./l_hog_gf(:)./H_target);
epsilon_dmax_mba_hog_gf=DR_mag_hog_gf      (:)./(1+H_target.*l_hog_gf(:)     .^2./18./I_fot_1./EGratio);
epsilon_dmax_mba_hog_bl=DR_mag_hog_fr_i_el (:)./(1+H_target.*l_hog_gf(:)     .^2./18./I_fot_1./EGratio);
 
if num_storey==0
Epsilon_h=eps_avg_hog_fr_i_el(:);
Epsilon_bmax=epsilon_bmax_mba_hog_bl;
Epsilon_dmax=epsilon_dmax_mba_hog_bl;
epsilon_br_mba_hog_bl=Epsilon_bmax+Epsilon_h;
epsilon_dr_mba_hog_bl=Epsilon_h.*((1-ni_foot)./2)+(Epsilon_h.^2.*((1+ni_foot)./2).^2+Epsilon_dmax.^2).^0.5;   

Epsilon_h=eps_avg_sag_fr_i_el(:);
Epsilon_bmax=epsilon_bmax_mba_sag_bl;
Epsilon_dmax=epsilon_dmax_mba_sag_bl;
epsilon_br_mba_sag_bl=Epsilon_bmax+Epsilon_h;
epsilon_dr_mba_sag_bl=Epsilon_h.*((1-ni_foot)./2)+(Epsilon_h.^2.*((1+ni_foot)./2).^2+Epsilon_dmax.^2).^0.5;


Epsilon_h=eps_avg_hog_gf(:);
Epsilon_bmax=epsilon_bmax_mba_hog_gf;
Epsilon_dmax=epsilon_dmax_mba_hog_gf;
epsilon_br_mba_hog_gf=Epsilon_bmax+Epsilon_h;
epsilon_dr_mba_hog_gf=Epsilon_h.*((1-ni_foot)./2)+(Epsilon_h.^2.*((1+ni_foot)./2).^2+Epsilon_dmax.^2).^0.5;   

Epsilon_h=eps_avg_sag_gf(:);
Epsilon_bmax=epsilon_bmax_mba_sag_gf;
Epsilon_dmax=epsilon_dmax_mba_sag_gf;
epsilon_br_mba_sag_gf=Epsilon_bmax+Epsilon_h;
epsilon_dr_mba_sag_gf=Epsilon_h.*((1-ni_foot)./2)+(Epsilon_h.^2.*((1+ni_foot)./2).^2+Epsilon_dmax.^2).^0.5;

end










Table_epsv1=[max(max(epsilon_br_top))   max(max(epsilon_br_bot))      max(max(epsilon_dr));
           max(epsilon_br_dba_hog_bl) max(epsilon_br_dba_sag_bl)    max(max(epsilon_dr_dba_sag_bl,epsilon_dr_dba_hog_bl));
           max(epsilon_br_mba_hog_bl) max(epsilon_br_mba_sag_bl)    max(max(epsilon_dr_mba_sag_bl,epsilon_dr_mba_hog_bl))        ] ;

if switch_tun_adv==1
[Table_eps(1,1),Table_i_ysc(1,1)]=max(max(epsilon_br_top));
[Table_eps(1,2),Table_i_ysc(1,2)]=max(max(epsilon_br_bot));
[Table_eps(1,3),Table_i_ysc(1,3)]=max(max(epsilon_dr));
[Table_eps(2,1),Table_i_ysc(2,1)]=max(epsilon_br_dba_hog_bl);
[Table_eps(2,2),Table_i_ysc(2,2)]=max(epsilon_br_dba_sag_bl);
[Table_eps(2,3),Table_i_ysc(2,3)]= max(max(epsilon_dr_dba_sag_bl,epsilon_dr_dba_hog_bl));
[Table_eps(3,1),Table_i_ysc(3,1)]=max(epsilon_br_mba_hog_bl);
[Table_eps(3,2),Table_i_ysc(3,2)]=max(epsilon_br_mba_sag_bl);
[Table_eps(3,3),Table_i_ysc(3,3)]= max(max(epsilon_dr_mba_sag_bl,epsilon_dr_mba_hog_bl));
Table_ysc =ys_vect(Table_i_ysc );  

vect_Table_eps=Table_eps(:);
vect_Table_ysc=Table_ysc(:);

result_TSI_strain=[result_TSI,vect_Table_eps',vect_Table_ysc'];

else
    
result_TSI_strain=[result_TSI,Table_epsv1(:)'];  
    
end






%   
% Legend1 =legend( [h(1),h(21),h(22),h(3),h(4),h(5),h(6),h(7),h(8),h(9),h(10)...
%     ,h(23),h(24),h(25),h(26),h(27),h(28),h(29),h(210)],...
%     {'epsilon-dmax','epsilon-bmax-top','epsilon-bmax-bot'...
%     ,'epsilon-dmax-dba-sag-tsi','epsilon-bmax-dba-sag-tsi','epsilon-dmax-dba-sag-gf','epsilon-bmax-dba-sag-gf','epsilon-dmax-dba-hog-tsi','epsilon-bmax-dba-hog-tsi'...
%     ,'epsilon-dmax-dba-hog-gf','epsilon-bmax-dba-hog-gf'...
%     ,'epsilon-dmax-mba-sag-tsi','epsilon-bmax-mba-sag-tsi','epsilon-dmax-mba-sag-gf','epsilon-bmax-mba-sag-gf','epsilon-dmax-mba-hog-tsi','epsilon-bmax-mba-hog-tsi'...
%     ,'epsilon-dmax-mba-hog-gf','epsilon-bmax-mba-hog-gf'},...
%     'Location','southeastoutside');%,'PlotBoxAspectRatio',ratio_leg21'Position',pos_leg_1','PlotBoxAspectRatioMode','manual'
% legend('boxon')





end
end


%%

% Ks_Jac=1/12*Efoot/Es*(dfoot/B_X_foot)^3;
% Delta_norm=max(Uzpg_fr_i_off)/B_X_foot/(Fz_foot_centr/bfoot/1000);
% result_jacok=[Ks_Jac;Delta_norm];

%% Something about the frame on continous beam(check)
if nfoot==1&&num_storey>0
    Elements_number=size(Element_nodes,1);
    for ii=1:Elements_number
        
        
        %Selection of element's Dofs
        ENs=Element_nodes(ii,:);  %Nodes of this element
        EDofs=[[1 2 3]+3*(ENs(1)-1) [1 2 3]+3*(ENs(2)-1)]; %Degress of freedom associated to this element
        
        % U_qbeam_fixhead_foot(EDofs)
        Fin_frame=KG_el_qbeam_fixhead_foot(:,:,ii)*U_qbeam_fixhead_foot;
        
        Fin_el(ii,:)=Fin_frame(EDofs)';
        %FL_qbeam_fixhead_foot ,KG_el_qbeam_fixhead_foot,Element_nodes
        
    end
end




if switch_plot==0
else

%% PLOT

ind=num_incr_TUNNEL;
figure(23);hold on;clf;
subplot(2,1,1); 
hold on;%xxso
h(1)=plot(xnod_tot,vect_Ugf_z(:,ind)*1000,'g-');

h(3)=plot(xnod_tot,vect_Utotf_fr_i_ep_deltaT_z(:,ind)*1000,'ro','linewidth',1,'MarkerSize',5);
h(2)=plot(xnod_tot,vect_Utotf_fr_i_el_deltaT_z(:,ind)*1000,'ko','linewidth',1,'MarkerSize',5);
h(4)=plot(xnod_tot,Uzpg_fr_i_off_deltaT(:,ind)*1000,'bx','linewidth',1,'MarkerSize',5);
h(5)=plot(xnod_tot,Uzpg_fr_i_deltaT(:,ind)    *1000,'k-x','linewidth',1,'MarkerSize',5);


% % % % % % % % % % % load inp_vlt_cd13_id09.mat
% % % % % % % % % % % load inp_new_mesh_all_cd13_id09.mat
% % % % % % % % % % % plot(75*X(1,:)/1000,75*DY{10}(1,:),'ko')



Legend1 =legend( [h(1),h(2),h(3),h(4),h(5)],...
    {'GF','GL-EL','GL-EP','NA-EL','GL-EL'},...
    'Location','southeastoutside');

ylabel('u_z (mm)');% 
xlabel('x_r (m)');% 
% hLegend =legend(h([10,1,3,5]) ,    {'STR (EL)','STR (EP)','GF','SLIDER'},'location','eastoutside');
set(gca,'FontName','Times New Roman','linewidth',1.0,'YDir','reverse','ticklength',1.0*get(gca,'ticklength'));%,'XLim',[-10 10],'YLim',[0 13],,'XTickLabel',XTickLabel_vector_empty,'yTick',yTick_vector'xTick',xTick_vector_z,'yTick',yTick_vector_z,'ticklength',2.0*get(gca,'ticklength'),,'XTickLabel',XTickLabel_vector_empty,,'XTickLabel',XTickLabel_vector_empty
title(['Tunnelling induced movements for Vlt: ' num2str(Vltp_vect(ind)) ' %. in (x_r,y_r,z_r)'])
box('on')

% % 
subplot(2,1,2);hold on;
% 
% h(20)=plot(xnod_tot/N,vect_Utotf_fr_i_el_P_x(:,ind)/N*1000,'b-','linewidth',2,'MarkerSize',10);

h(1)=plot(xnod_tot,vect_Ugf_x(:,ind)*1000,'g-');
h(2)=plot(xnod_tot,vect_Utotf_fr_i_el_deltaT_x(:,ind)*1000,'kx','linewidth',1,'MarkerSize',5);
h(3)=plot(xnod_tot,vect_Utotf_fr_i_ep_deltaT_x(:,ind)*1000,'rx','linewidth',1,'MarkerSize',5);
h(4)=plot(xnod_tot,Uxpg_fr_i_off_deltaT(:,ind)*1000,'bo','linewidth',1,'MarkerSize',5);
h(5)=plot(xnod_tot,Uxpg_fr_i_deltaT(:,ind)    *1000,'k-o','linewidth',1,'MarkerSize',5);
% 
% h(4)=plot(xnod_tot,Uxpg_fr_i_off(:)*1000,'b:o','linewidth',1,'MarkerSize',5);
% h(5)=plot(xnod_tot,Uxpg_fr_i(:)    *1000,'k:s','linewidth',1,'MarkerSize',5);

% % % % % % % % plot(75*X(1,:)/1000,75*DX{10}(1,:),'ko')

Legend1 =legend( [h(1),h(2),h(3),h(4),h(5)],...
    {'GF','GL-EL','GL-EP','NA-EL','GL-EL'},...
    'Location','southeastoutside');

ylabel('u_x (mm)');% 
xlabel('x_r (m)');% 
% % hLegend =legend(h([20,2,4,6]) ,    {'STR (EL)','STR (EP)','GF','SLIDER'},'location','eastoutside');
set(gca,'FontName','Times New Roman','linewidth',1.0,'YDir','reverse','ticklength',1.0*get(gca,'ticklength'));%,'XLim',[-10 10],'YLim',[0 2],,'XTickLabel',XTickLabel_vector_empty,'yTick',yTick_vector'xTick',xTick_vector_z,'yTick',yTick_vector_z,'ticklength',2.0*get(gca,'ticklength'),,'XTickLabel',XTickLabel_vector_empty,,'XTickLabel',XTickLabel_vector_empty
box('on')
% % grid on
% 

saveas(gcf,'Vlt50.png')


ind=num_incr_TUNNEL/5;
figure(230);hold on;clf;
subplot(2,1,1); 
hold on;%xxso
h(1)=plot(xnod_tot,vect_Ugf_z(:,ind)*1000,'g-');

h(3)=plot(xnod_tot,vect_Utotf_fr_i_ep_deltaT_z(:,ind)*1000,'ro','linewidth',1,'MarkerSize',5);
h(2)=plot(xnod_tot,vect_Utotf_fr_i_el_deltaT_z(:,ind)*1000,'ko','linewidth',1,'MarkerSize',5);
h(4)=plot(xnod_tot,Uzpg_fr_i_off_deltaT(:,ind)*1000,'bx','linewidth',1,'MarkerSize',5);
h(5)=plot(xnod_tot,Uzpg_fr_i_deltaT(:,ind)    *1000,'k-x','linewidth',1,'MarkerSize',5);


% % % % % % % % % % % load inp_vlt_cd13_id09.mat
% % % % % % % % % % % load inp_new_mesh_all_cd13_id09.mat
% % % % % % % % % % % plot(75*X(1,:)/1000,75*DY{10}(1,:),'ko')



Legend1 =legend( [h(1),h(2),h(3),h(4),h(5)],...
    {'GF','GL-EL','GL-EP','NA-EL','GL-EL'},...
    'Location','southeastoutside');

ylabel('u_z (mm)');% 
xlabel('x_r (m)');% 
% hLegend =legend(h([10,1,3,5]) ,    {'STR (EL)','STR (EP)','GF','SLIDER'},'location','eastoutside');
set(gca,'FontName','Times New Roman','linewidth',1.0,'YDir','reverse','ticklength',1.0*get(gca,'ticklength'));%,'XLim',[-10 10],'YLim',[0 13],,'XTickLabel',XTickLabel_vector_empty,'yTick',yTick_vector'xTick',xTick_vector_z,'yTick',yTick_vector_z,'ticklength',2.0*get(gca,'ticklength'),,'XTickLabel',XTickLabel_vector_empty,,'XTickLabel',XTickLabel_vector_empty
title(['Tunnelling induced movements for Vlt: ' num2str(Vltp_vect(ind)) ' %. in (x_r,y_r,z_r)'])
box('on')

% % 
subplot(2,1,2);hold on;
% 
% h(20)=plot(xnod_tot/N,vect_Utotf_fr_i_el_P_x(:,ind)/N*1000,'b-','linewidth',2,'MarkerSize',10);

h(1)=plot(xnod_tot,vect_Ugf_x(:,ind)*1000,'g-');
h(2)=plot(xnod_tot,vect_Utotf_fr_i_el_deltaT_x(:,ind)*1000,'kx','linewidth',1,'MarkerSize',5);
h(3)=plot(xnod_tot,vect_Utotf_fr_i_ep_deltaT_x(:,ind)*1000,'rx','linewidth',1,'MarkerSize',5);
h(4)=plot(xnod_tot,Uxpg_fr_i_off_deltaT(:,ind)*1000,'bo','linewidth',1,'MarkerSize',5);
h(5)=plot(xnod_tot,Uxpg_fr_i_deltaT(:,ind)    *1000,'k-o','linewidth',1,'MarkerSize',5);
% 
% h(4)=plot(xnod_tot,Uxpg_fr_i_off(:)*1000,'b:o','linewidth',1,'MarkerSize',5);
% h(5)=plot(xnod_tot,Uxpg_fr_i(:)    *1000,'k:s','linewidth',1,'MarkerSize',5);

% % % % % % % % plot(75*X(1,:)/1000,75*DX{10}(1,:),'ko')

Legend1 =legend( [h(1),h(2),h(3),h(4),h(5)],...
    {'GF','GL-EL','GL-EP','NA-EL','GL-EL'},...
    'Location','southeastoutside');

ylabel('u_x (mm)');% 
xlabel('x_r (m)');% 
% % hLegend =legend(h([20,2,4,6]) ,    {'STR (EL)','STR (EP)','GF','SLIDER'},'location','eastoutside');
set(gca,'FontName','Times New Roman','linewidth',1.0,'YDir','reverse','ticklength',1.0*get(gca,'ticklength'));%,'XLim',[-10 10],'YLim',[0 2],,'XTickLabel',XTickLabel_vector_empty,'yTick',yTick_vector'xTick',xTick_vector_z,'yTick',yTick_vector_z,'ticklength',2.0*get(gca,'ticklength'),,'XTickLabel',XTickLabel_vector_empty,,'XTickLabel',XTickLabel_vector_empty
box('on')
% % grid on
% 

saveas(gcf,'Vlt10.png')

ind=num_incr_TUNNEL/2;
figure(2300);hold on;clf;
subplot(2,1,1); 
hold on;%xxso
h(1)=plot(xnod_tot,vect_Ugf_z(:,ind)*1000,'g-');

h(3)=plot(xnod_tot,vect_Utotf_fr_i_ep_deltaT_z(:,ind)*1000,'ro','linewidth',1,'MarkerSize',5);
h(2)=plot(xnod_tot,vect_Utotf_fr_i_el_deltaT_z(:,ind)*1000,'ko','linewidth',1,'MarkerSize',5);
h(4)=plot(xnod_tot,Uzpg_fr_i_off_deltaT(:,ind)*1000,'bx','linewidth',1,'MarkerSize',5);
h(5)=plot(xnod_tot,Uzpg_fr_i_deltaT(:,ind)    *1000,'k-x','linewidth',1,'MarkerSize',5);


% % % % % % % % % % % load inp_vlt_cd13_id09.mat
% % % % % % % % % % % load inp_new_mesh_all_cd13_id09.mat
% % % % % % % % % % % plot(75*X(1,:)/1000,75*DY{10}(1,:),'ko')



Legend1 =legend( [h(1),h(2),h(3),h(4),h(5)],...
    {'GF','GL-EL','GL-EP','NA-EL','GL-EL'},...
    'Location','southeastoutside');

ylabel('u_z (mm)');% 
xlabel('x_r (m)');% 
% hLegend =legend(h([10,1,3,5]) ,    {'STR (EL)','STR (EP)','GF','SLIDER'},'location','eastoutside');
set(gca,'FontName','Times New Roman','linewidth',1.0,'YDir','reverse','ticklength',1.0*get(gca,'ticklength'));%,'XLim',[-10 10],'YLim',[0 13],,'XTickLabel',XTickLabel_vector_empty,'yTick',yTick_vector'xTick',xTick_vector_z,'yTick',yTick_vector_z,'ticklength',2.0*get(gca,'ticklength'),,'XTickLabel',XTickLabel_vector_empty,,'XTickLabel',XTickLabel_vector_empty
title(['Tunnelling induced movements for Vlt: ' num2str(Vltp_vect(ind)) ' %. in (x_r,y_r,z_r)'])
box('on')

% % 
subplot(2,1,2);hold on;
% 
% h(20)=plot(xnod_tot/N,vect_Utotf_fr_i_el_P_x(:,ind)/N*1000,'b-','linewidth',2,'MarkerSize',10);

h(1)=plot(xnod_tot,vect_Ugf_x(:,ind)*1000,'g-');
h(2)=plot(xnod_tot,vect_Utotf_fr_i_el_deltaT_x(:,ind)*1000,'kx','linewidth',1,'MarkerSize',5);
h(3)=plot(xnod_tot,vect_Utotf_fr_i_ep_deltaT_x(:,ind)*1000,'rx','linewidth',1,'MarkerSize',5);
h(4)=plot(xnod_tot,Uxpg_fr_i_off_deltaT(:,ind)*1000,'bo','linewidth',1,'MarkerSize',5);
h(5)=plot(xnod_tot,Uxpg_fr_i_deltaT(:,ind)    *1000,'k-o','linewidth',1,'MarkerSize',5);
% 
% h(4)=plot(xnod_tot,Uxpg_fr_i_off(:)*1000,'b:o','linewidth',1,'MarkerSize',5);
% h(5)=plot(xnod_tot,Uxpg_fr_i(:)    *1000,'k:s','linewidth',1,'MarkerSize',5);

% % % % % % % % plot(75*X(1,:)/1000,75*DX{10}(1,:),'ko')

Legend1 =legend( [h(1),h(2),h(3),h(4),h(5)],...
    {'GF','GL-EL','GL-EP','NA-EL','GL-EL'},...
    'Location','southeastoutside');

ylabel('u_x (mm)');% 
xlabel('x_r (m)');% 
% % hLegend =legend(h([20,2,4,6]) ,    {'STR (EL)','STR (EP)','GF','SLIDER'},'location','eastoutside');
set(gca,'FontName','Times New Roman','linewidth',1.0,'YDir','reverse','ticklength',1.0*get(gca,'ticklength'));%,'XLim',[-10 10],'YLim',[0 2],,'XTickLabel',XTickLabel_vector_empty,'yTick',yTick_vector'xTick',xTick_vector_z,'yTick',yTick_vector_z,'ticklength',2.0*get(gca,'ticklength'),,'XTickLabel',XTickLabel_vector_empty,,'XTickLabel',XTickLabel_vector_empty
box('on')
% % grid on
% 
saveas(gcf,'Vlt25.png')


%%%%%%%%%%%%%
%PLOT 3D FIGURES IN THE CASE OF TUNNEL ADVANCEMENTS
%Plot results for the latest tunnel advancement
if switch_tun_adv==1
    % Ground grid
    dgrid=1;
    flagclay=0;
    meshinput=-50:dgrid:50;
    [Xgr,Ygr]=meshgrid(meshinput);
    xgr=Xgr(:);
    ygr=Ygr(:);
    
    
    % Building grid auxillary
    [dxbg2,dybg2,dzbg2]=comp3_AF(xbg-dgrid,ybg,0,ys,yf,y0,Kx,Ky,ht,Vltp/100,Rt*2);
    [dxbg3,dybg3,dzbg3]=comp3_AF(xbg,ybg-dgrid,0,ys,yf,y0,Kx,Ky,ht,Vltp/100,Rt*2);
    
    % For grid
    [dxgr,dygr,dzgr]=comp3_AF(xgr,ygr,0,ys,yf,y0,Kx,Ky,ht,Vltp/100,Rt*2);
    
    if switch_tun_adv_twin==1
    % For grid
    
%         xgr_twin=+(forig_twin*cos(theta_twin-theta)-forig)          +xgr*cos(theta-theta_twin)+ygr*sin(theta-theta_twin);
%         ygr_twin=-forig_twin*sin(theta_twin-theta)                  -xgr*sin(theta-theta_twin)+ygr*cos(theta-theta_twin);
        xgr_twin=(xgr+(forig_twin*cos(theta_twin-theta)-forig))          ;
        ygr_twin=(ygr-forig_twin*sin(theta_twin-theta) )                 ;
        xgr_twin_rot=         +xgr_twin*cos(theta-theta_twin)+ygr_twin*sin(theta-theta_twin);
        ygr_twin_rot=         -xgr_twin*sin(theta-theta_twin)+ygr_twin*cos(theta-theta_twin);


        
        
        [dxgr_twin,dygr_twin,dzgr_twin]=comp3_AF(xgr_twin_rot,ygr_twin_rot,0,ys_twin,yf_twin,y0_twin,Kx_twin,Ky_twin,ht_twin,Vltp_twin/100,Rt_twin*2);


        %I need to project the displacements from the twin refer system to
        %the first tunnel reference system
        dxgr_twin_proj=+dxgr_twin*cos(theta-theta_twin)-dygr_twin*sin(theta-theta_twin);
        dygr_twin_proj=+dxgr_twin*sin(theta-theta_twin)+dygr_twin*cos(theta-theta_twin);
        dzgr_twin_proj=dzgr_twin;

    
        dxgr=[dxgr+dxgr_twin_proj];
        dygr=[dygr+dygr_twin_proj];
        dzgr=[dzgr+dzgr_twin_proj];
    end
    
    
    
    % Calc strain
    epsxbl=diff(dxbl)/(xbl(2)-xbl(1));
    epsxblp=(dxbg-dxbg2)/dgrid;
    epsyblp=(dybg-dybg3)/dgrid;
    epsxyblp=(dxbg-dxbg3)/dgrid/2+...
        (dybg-dybg2)/dgrid/2;
    xbg = xbg(:);
    ybg = ybg(:);
    dxbg= dxbg(:);
    dybg= dybg(:);
    dzbg= dzbg(:);
   
    
    dxbg_SSI=+cos(theta)*vect_Utotf_fr_i_el_deltaT_x(:,end)-dybl*sin(theta);%Need to go from local to global
    dybg_SSI=+sin(theta)*vect_Utotf_fr_i_el_deltaT_x(:,end)+dybl*cos(theta);%Need to go from local to global
    dzbg_SSI=vect_Utotf_fr_i_el_deltaT_z(:,end);
    
 
    
    
    figure
    subplot(2,2,1)
    scatter3(xgr,ygr,dxgr*1e3,10,'bx')
    hold on
    scatter3(xbg,ybg,dxbg    *1e3,50,'ro','filled')
    scatter3(xbg,ybg,dxbg_SSI*1e3,50,'co','filled')
    title('\deltax GF displacements (mm) - Global reference system (x_g,y_g,z_g)')
    xlabel('x'); ylabel('y');
        subplot(2,2,2)
    scatter3(xgr,ygr,dygr*1e3,10,'bx')
    hold on
    scatter3(xbg,ybg,dybg*1e3,50,'ro','filled')
    scatter3(xbg,ybg,dybg_SSI*1e3,50,'co','filled')
    title('\deltay GF displacements (mm) - Global reference system (x_g,y_g,z_g)')
    xlabel('x'); ylabel('y');
    subplot(2,2,3)
    scatter3(xgr,ygr,dzgr*1e3,10,'bx')
    hold on
    scatter3(xbg,ybg,dzbg    *1e3,50,'ro','filled')
    scatter3(xbg,ybg,dzbg_SSI*1e3,50,'co','filled')
    title('\deltaz GF displacements (mm) - Global reference system (x_g,y_g,z_g) ')
    xlabel('x'); ylabel('y');
    subplot(2,2,4)
    scatter(xgr,ygr,10,'bx')
    hold on
    scatter(zeros(length(ys:dgrid:yf),1),ys:dgrid:yf,50,'go','filled')
    scatter(xbg,ybg,50,'ro','filled')
    title('Building position and tunnel advance - Global reference system (x_g,y_g,z_g)')
    legend('Ground','Tunnel','Building')
    xlabel('x'); ylabel('y');
    
    figure
    subplot(1,2,1)
    surfc(Xgr,Ygr,reshape(dzgr,[length(meshinput),length(meshinput)]),...
        'EdgeColor','none')
    title('3D GF settlement trough')
    ylabel('Y (m)')
    xlabel('X (m)')
    zlabel('\DeltaZ (mm)')
    subplot(1,2,2)
    scatter(xgr,ygr,10,'bx')
    hold on
    scatter(zeros(length(ys:dgrid:yf),1),ys:dgrid:yf,50,'go','filled')
    scatter(xbg,ybg,50,'ro','filled')
    ylabel('Y (m)')
    xlabel('X (m)')
    title('Building position and tunnel advance')
    
%     figure
%     subplot(2,2,1)
%     title('Camos Figure 8 strain validation')
%     plot(abs(xbl(2:end)-max(xbl(2:end))),epsxbl)
%     ylabel('\epsilon_h (mm)')
%     xlabel('r (m)')
%     subplot(2,2,2)
%     plot(abs(xbl-max(xbl)),epsxblp)
%     ylabel('\epsilon_{h,xx}  (mm)')
%     xlabel('r (m)')
%     subplot(2,2,3)
%     plot(abs(xbl-max(xbl)),epsyblp)
%     ylabel('\epsilon_{h,yy} (mm)')
%     xlabel('\r (m)')
%     subplot(2,2,4)
%     plot(abs(xbl-max(xbl)),epsxyblp)
%     ylabel('\epsilon_{h,xy} (mm)')
%     xlabel('\r (m)')
    
    figure
    subplot(2,2,1)
    plot(xbl,dxbl*1e3)
    ylabel('\Deltax" (mm)')
    xlabel('x" (m)')
    subplot(2,2,2)
    plot(xbl,dybl*1e3)
    ylabel('\Deltay" (mm)')
    xlabel('x" (m)')
    subplot(2,2,3)
    plot(xbl,dzbl*1e3)
    ylabel('\Deltaz" (mm)')
    xlabel('x" (m)')
    subplot(2,2,4)
    scatter(xgr,ygr,10,'b.')
    hold on
    scatter(zeros(length(ys:dgrid:yf),1),ys:dgrid:yf,50,'go','filled')
    scatter(xbg,ybg,50,'ro','filled')
    title('Building position and tunnel advance')
    legend('Ground','Tunnel','Building')
end


if switch_tun_adv==1
figure(31);hold on;
title('epsb and epsd')
h(1)=plot(ys_vect,max(epsilon_br_top),'-r');
h(2)=plot(ys_vect,max(epsilon_br_bot),':b');
h(3)=plot(ys_vect,max(epsilon_dr),'-c');

h(10)=plot(ys_vect,(epsilon_br_dba_hog_bl),'or');
h(20)=plot(ys_vect,(epsilon_br_dba_sag_bl),'ob');
h(30)=plot(ys_vect,max(epsilon_dr_dba_sag_bl,epsilon_dr_dba_hog_bl),'oc');

h(10)=plot(ys_vect,(epsilon_br_mba_hog_bl),'xr');
h(20)=plot(ys_vect,(epsilon_br_mba_sag_bl),'xb');
h(30)=plot(ys_vect,(max(epsilon_dr_mba_sag_bl,epsilon_dr_mba_hog_bl)),'xc');

else
    figure(31);hold on;
title('epsb and epsd')
h(1)=plot(Vltp_vect,max(epsilon_br_top),'-b');
h(2)=plot(Vltp_vect,max(epsilon_br_bot),':b');
h(3)=plot(Vltp_vect,max(epsilon_dr),'-c');

h(10)=plot(Vltp_vect,(epsilon_br_dba_hog_bl),'-ob');
h(20)=plot(Vltp_vect,(epsilon_br_dba_sag_bl),':ob');
h(30)=plot(Vltp_vect,max(epsilon_dr_dba_sag_bl,epsilon_dr_dba_hog_bl),'-oc');

h(10)=plot(Vltp_vect,(epsilon_br_mba_hog_bl),'-xb');
h(20)=plot(Vltp_vect,(epsilon_br_mba_sag_bl),':xb');
h(30)=plot(Vltp_vect,(max(epsilon_dr_mba_sag_bl,epsilon_dr_mba_hog_bl)),'-xc');

end



figure(21);hold on;
title('epsb and epsd')
h(1)=plot(Xord,epsilon_dmax(:,end),'-k');
h(21)=plot(Xord,epsilon_bmax_top(:,end),'-r');
h(22)=plot(Xord,epsilon_bmax_bot(:,end),'-m');

h(3)=plot([Xord(1),Xord(end)],[epsilon_dmax_bur_sag_bl(end),epsilon_dmax_bur_sag_bl(end)],'-xk');
h(4)=plot([Xord(1),Xord(end)],[epsilon_bmax_bur_sag_bl(end),epsilon_bmax_bur_sag_bl(end)],'-xm');
 
h(5)=plot([Xord(1),Xord(end)],[epsilon_dmax_bur_sag_gf(end),epsilon_dmax_bur_sag_gf(end)],'--xk');
h(6)=plot([Xord(1),Xord(end)],[epsilon_bmax_bur_sag_gf(end),epsilon_bmax_bur_sag_gf(end)],'--xm');

h(7)=plot([Xord(1),Xord(end)],[epsilon_dmax_bur_hog_bl(end),epsilon_dmax_bur_hog_bl(end)],'-ok');
h(8)=plot([Xord(1),Xord(end)],[epsilon_bmax_bur_hog_bl(end),epsilon_bmax_bur_hog_bl(end)],'-or');

h(9)=plot([Xord(1),Xord(end)],[epsilon_dmax_bur_hog_gf(end),epsilon_dmax_bur_hog_gf(end)],'--ok');
h(10)=plot([Xord(1),Xord(end)],[epsilon_bmax_bur_hog_gf(end),epsilon_bmax_bur_hog_gf(end)],'--or');



h(23)=plot([Xord(1),Xord(end)],[epsilon_dmax_mba_sag_bl(end),epsilon_dmax_mba_sag_bl(end)],'-+k','MarkerSize',15);
h(24)=plot([Xord(1),Xord(end)],[epsilon_bmax_mba_sag_bl(end),epsilon_bmax_mba_sag_bl(end)],'-+m','MarkerSize',15);
 
h(25)=plot([Xord(1),Xord(end)],[epsilon_dmax_mba_sag_gf(end),epsilon_dmax_mba_sag_gf(end)],'--+k','MarkerSize',15);
h(26)=plot([Xord(1),Xord(end)],[epsilon_bmax_mba_sag_gf(end),epsilon_bmax_mba_sag_gf(end)],'--+m','MarkerSize',15);

h(27)=plot([Xord(1),Xord(end)],[epsilon_dmax_mba_hog_bl(end),epsilon_dmax_mba_hog_bl(end)],'-sk','MarkerSize',15);
h(28)=plot([Xord(1),Xord(end)],[epsilon_bmax_mba_hog_bl(end),epsilon_bmax_mba_hog_bl(end)],'-sr','MarkerSize',15);

h(29)=plot([Xord(1),Xord(end)],[epsilon_dmax_mba_hog_gf(end),epsilon_dmax_mba_hog_gf(end)],'--sk','MarkerSize',15);
h(210)=plot([Xord(1),Xord(end)],[epsilon_bmax_mba_hog_gf(end),epsilon_bmax_mba_hog_gf(end)],'--sr','MarkerSize',15);


  
Legend1 =legend( [h(1),h(21),h(22),h(3),h(4),h(5),h(6),h(7),h(8),h(9),h(10)...
    ,h(23),h(24),h(25),h(26),h(27),h(28),h(29),h(210)],...
    {'epsilon-dmax','epsilon-bmax-top','epsilon-bmax-bot'...
    ,'epsilon-dmax-dba-sag-tsi','epsilon-bmax-dba-sag-tsi','epsilon-dmax-dba-sag-gf','epsilon-bmax-dba-sag-gf','epsilon-dmax-dba-hog-tsi','epsilon-bmax-dba-hog-tsi'...
    ,'epsilon-dmax-dba-hog-gf','epsilon-bmax-dba-hog-gf'...
    ,'epsilon-dmax-mba-sag-tsi','epsilon-bmax-mba-sag-tsi','epsilon-dmax-mba-sag-gf','epsilon-bmax-mba-sag-gf','epsilon-dmax-mba-hog-tsi','epsilon-bmax-mba-hog-tsi'...
    ,'epsilon-dmax-mba-hog-gf','epsilon-bmax-mba-hog-gf'},...
    'Location','southeastoutside');%,'PlotBoxAspectRatio',ratio_leg21'Position',pos_leg_1','PlotBoxAspectRatioMode','manual'
legend('boxon')

figure(210);hold on;
title('epsbt and epsdt')
% h(1)=plot(Xord,epsilon_br,'-b');
h(11)=plot(Xord,epsilon_br_top(:,end),'-b');
h(12)=plot(Xord,epsilon_br_bot(:,end),':b');
h(2)=plot(Xord,epsilon_dr(:,end),'-c');
h(3)=plot([Xord(1),Xord(end)],[epsilon_dr_dba_sag_bl(end),epsilon_dr_dba_sag_bl(end)],'-xc');
h(4)=plot([Xord(1),Xord(end)],[epsilon_br_dba_sag_bl(end),epsilon_br_dba_sag_bl(end)],'-xb');

h(5)=plot([Xord(1),Xord(end)],[epsilon_dr_dba_hog_bl(end),epsilon_dr_dba_hog_bl(end)],'-oc');
h(6)=plot([Xord(1),Xord(end)],[epsilon_br_dba_hog_bl(end),epsilon_br_dba_hog_bl(end)],'-ob');

% h(43)=plot([Xord(1),Xord(end)],[epsilon_dr_dba_sag_gf(end),epsilon_dr_dba_sag_gf(end)],'--xc');
% h(44)=plot([Xord(1),Xord(end)],[epsilon_br_dba_sag_gf(end),epsilon_br_dba_sag_gf(end)],'--xb');
% 
% h(45)=plot([Xord(1),Xord(end)],[epsilon_dr_dba_hog_gf(end),epsilon_dr_dba_hog_gf(end)],'--oc');
% h(46)=plot([Xord(1),Xord(end)],[epsilon_br_dba_hog_gf(end),epsilon_br_dba_hog_gf(end)],'--ob');




h(23)=plot([Xord(1),Xord(end)],[epsilon_dr_mba_sag_bl(end),epsilon_dr_mba_sag_bl(end)],'-+c','MarkerSize',15);
h(24)=plot([Xord(1),Xord(end)],[epsilon_br_mba_sag_bl(end),epsilon_br_mba_sag_bl(end)],'-+b','MarkerSize',15);

h(25)=plot([Xord(1),Xord(end)],[epsilon_dr_mba_hog_bl(end),epsilon_dr_mba_hog_bl(end)],'-sc','MarkerSize',15);
h(26)=plot([Xord(1),Xord(end)],[epsilon_br_mba_hog_bl(end),epsilon_br_mba_hog_bl(end)],'-sb','MarkerSize',15);

% h(423)=plot([Xord(1),Xord(end)],[epsilon_dr_mba_sag_gf(end),epsilon_dr_mba_sag_gf(end)],'--+c','MarkerSize',15);
% h(424)=plot([Xord(1),Xord(end)],[epsilon_br_mba_sag_gf(end),epsilon_br_mba_sag_gf(end)],'--+b','MarkerSize',15);
% 
% h(425)=plot([Xord(1),Xord(end)],[epsilon_dr_mba_hog_gf(end),epsilon_dr_mba_hog_gf(end)],'--sc','MarkerSize',15);
% h(426)=plot([Xord(1),Xord(end)],[epsilon_br_mba_hog_gf(end),epsilon_br_mba_hog_gf(end)],'--sb','MarkerSize',15);
% 


Legend1 =legend( [h(11),h(12),h(2),h(3),h(4),h(5),h(6),h(23),h(24),h(25),h(26)],...
    {'epsilon-bt-top','epsilon-bt-bot','epsilon-dt',...
    'epsilon-dt-dba-sag-tsi','epsilon-bt-dba-sag-tsi','epsilon-dt-dba-hog-tsi','epsilon-bt-dba-hog-tsi',...
    'epsilon-dt-mba-sag-tsi','epsilon-bt-mba-sag-tsi','epsilon-dt-mba-hog-tsi','epsilon-bt-mba-hog-tsi'},...
    'Location','southeastoutside');%,'PlotBoxAspectRatio',ratio_leg21'Position',pos_leg_1','PlotBoxAspectRatioMode','manual'
legend('boxon')




end




