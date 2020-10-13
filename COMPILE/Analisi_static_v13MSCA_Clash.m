%%

%SCRIPT Analisi_static_v03MSCA_Clash
    %The CE solution is coupled to predict flim varying for a pile group
    %depending on vlt
        %Required action
        %Generalise the CE inputs with an excel spreadsheet
        %I need to modify flim in the horizontal direction
        %I need to check interaction factors
        %I need different limit force at tensile for shaft

%SCRIPT Analisi_static_v04MSCA_Clash
        %Started the implementation of a hyperbolic stress–strain law for
        %the local stiffness terms
        
%SCRIPT Analisi_static_v05MSCA_Clash  
        %Implemented with generical loading path with
        %switch_displfoundation

%SCRIPT Analisi_static_v08MSCA_Clash  
        %Implemented multilayered soil

%SCRIPT Analisi_static_v09MSCA_Clash  
        %Added sigmax sigmay and horizontal sliders
        
%SCRIPT Analisi_static_v10MSCA_Clash  
        %Added switch_soil=3 for directly input the Greenfield movement from measurement
%% ELABORATION INPUT DATA
% DESCRIPTIVE TEXT

%Position
x_centre=max(X)/2+min(X)/2;

%DISCRETIZZAZIONE DEL PALO
z=0:h_el:L;

% %%%%%%%%%%%%%
% %CHEN
z(1)=0.0001; 
ES=Es;
Es=Es(1);
ni=nis;
nis=nis(1);
% %%%%%%%%%%%%%



spessori=diff(z);
nnodes=length(z);
nstrati=length(spessori);       %numero degli strati
dz=spessori;                    %vettore della lunghezza degli elementi del singolo strato
npali=length(X);                %numero dei pali
dimK=nnodes*6;                           %dimensione matrice singolo palo
nelementi=L/h_el;                    %numero degli elementi del singolo palo



for i=1:nnodes;
	for ip=1:npali;
       z_global(i,ip)=z(i);
       x_global(i,ip)=X(ip);
       y_global(i,ip)=Y(ip);
	end
end

% %SPATIAL COORDINATES
% for i=1:length(zgauss);
% 	for ip=1:npali;
%        z_gauss_global(i,ip)=zgauss(i);
%        x_gauss_global(i,ip)=X(ip);
% 	end
% end



for i=1:npali;
    znod_tot(1+(nelementi+1)*(i-1):(nelementi+1)+(nelementi+1)*(i-1),1)=z';
    xnod_tot(1+(nelementi+1)*(i-1):1:(nelementi+1)+(nelementi+1)*(i-1),1)=X(i)+z*tan(pi./180*beta(1,i)) ;
    ynod_tot(1+(nelementi+1)*(i-1):1:(nelementi+1)+(nelementi+1)*(i-1),1)=Y(i)+z*tan(pi./180*alfa(1,i));
   
    Es_nod(1+(nelementi+1)*(i-1):(nelementi+1)+(nelementi+1)*(i-1),1)=Es(:,1);
    ni_nod(1+(nelementi+1)*(i-1):1:(nelementi+1)+(nelementi+1)*(i-1),1)=nis;
    
    Eb_nod_base(i,1)=Eb;
    nib_nod_base(i,1)=nib;
    
end
Gs_nod=Es_nod/2/(1+nis);
Gb_nod_base=Eb_nod_base/2/(1+nib);
Zi_nod=znod_tot;
Xi_nod=xnod_tot;
Yi_nod=ynod_tot;


for jb=1:npali
    Zi_nod_base(jb,1)=Zi_nod(jb*(nelementi+1));
    Xi_nod_base(jb,1)=Xi_nod(jb*(nelementi+1));
    Yi_nod_base(jb,1)=Yi_nod(jb*(nelementi+1));
end



%% CONDENSED STIFFNESS MATRIX of the STRUCTURE
% DESCRIPTIVE TEXT

if num_bay>0
    % %%CONDENSED MATRIX OF BUILODING ON FIXED CONSTRAINTS
    [ RR_fixhead ] = f_frame2D_cond_tie_pile_head_v1( num_bay,num_storey,spanx_frame,spanz_frame,Eframe,bc,dc,bb,db,bf,df,foundn );

    % %%CONDENSED MATRIX OF BUILODING ON HINGED CONSTRAINTS
    [ RR_hinhead ] = f_frame2D_cond_hinge_pile_head_v1( num_bay,num_storey,spanx_frame,spanz_frame,Eframe,bc,dc,bb,db,bf,df,foundn );

    % %%% EQUIVALENT BEAM
    % % DESCRIPTIVE TEXT
    % 
    [ EI_mat,EA_mat ] = f_beam_eq_coeff( ht,X,X(1),X(length(X)),Eframe,bf,df,bc,dc,bb,db,spanz_frame,spanx_frame,foundn,num_storey );
    % 
    % %tie_pile_head
    [ RR_fixhead_eq_beam] = f_beam2D_cond_tie_pile_head_varEI_v1( X,zeros(1,npali),EI_mat,EA_mat );
    % 
end
% The stiffness matrix computed by f_frame2D function has for column the
% reaction forces at pile heads due to unit constrained displacement. This
% condensed stiffness matrix need to be reorganized
KKframe_fixed=zeros(6*npali*nnodes);
KKframe_hinged=zeros(6*npali*nnodes);
if num_bay>0
    for ip=1:npali

        KKframe_fixed(1:(nelementi+1)*6:npali*(nelementi+1)*6,(ip-1)*6*(nelementi+1)+1)=RR_fixhead(1:3:npali*3,(ip-1)*3+1);
        KKframe_fixed(3:(nelementi+1)*6:npali*(nelementi+1)*6,(ip-1)*6*(nelementi+1)+1)=RR_fixhead(2:3:npali*3,(ip-1)*3+1);
        KKframe_fixed(5:(nelementi+1)*6:npali*(nelementi+1)*6,(ip-1)*6*(nelementi+1)+1)=RR_fixhead(3:3:npali*3,(ip-1)*3+1);

        KKframe_fixed(1:(nelementi+1)*6:npali*(nelementi+1)*6,(ip-1)*6*(nelementi+1)+3)=RR_fixhead(1:3:npali*3,(ip-1)*3+2);
        KKframe_fixed(3:(nelementi+1)*6:npali*(nelementi+1)*6,(ip-1)*6*(nelementi+1)+3)=RR_fixhead(2:3:npali*3,(ip-1)*3+2);
        KKframe_fixed(5:(nelementi+1)*6:npali*(nelementi+1)*6,(ip-1)*6*(nelementi+1)+3)=RR_fixhead(3:3:npali*3,(ip-1)*3+2);

        KKframe_fixed(1:(nelementi+1)*6:npali*(nelementi+1)*6,(ip-1)*6*(nelementi+1)+5)=RR_fixhead(1:3:npali*3,(ip-1)*3+3);
        KKframe_fixed(3:(nelementi+1)*6:npali*(nelementi+1)*6,(ip-1)*6*(nelementi+1)+5)=RR_fixhead(2:3:npali*3,(ip-1)*3+3);
        KKframe_fixed(5:(nelementi+1)*6:npali*(nelementi+1)*6,(ip-1)*6*(nelementi+1)+5)=RR_fixhead(3:3:npali*3,(ip-1)*3+3);
    end

    for ip=1:npali

        KKframe_hinged(1:(nelementi+1)*6:npali*(nelementi+1)*6,(ip-1)*6*(nelementi+1)+1)=RR_hinhead(1:2:npali*2,(ip-1)*2+1);
        KKframe_hinged(3:(nelementi+1)*6:npali*(nelementi+1)*6,(ip-1)*6*(nelementi+1)+1)=RR_hinhead(2:2:npali*2,(ip-1)*2+1);

        KKframe_hinged(1:(nelementi+1)*6:npali*(nelementi+1)*6,(ip-1)*6*(nelementi+1)+3)=RR_hinhead(1:2:npali*2,(ip-1)*2+2);
        KKframe_hinged(3:(nelementi+1)*6:npali*(nelementi+1)*6,(ip-1)*6*(nelementi+1)+3)=RR_hinhead(2:2:npali*2,(ip-1)*2+2);

    end
end




%% debug

% % KKframe_fixed=KKframe_hinged;

%% GREENFIELD DISPLACEMENT
% DESCRIPTIVE TEXT


if switch_soil==1
    [ Uffx,Uffz ]=u_LON(z_global,x_global,ht,Rt,epsilon,0,0.49999);
    Uffx=coeffux*Uffx;
    Uffz=coeffuz*Uffz;
    Uffy=zeros(nnodes,npali);
    
    [ Uffx_crown,Uffz_crown ]=u_LON(ht-Rt,0,ht,Rt,epsilon,0,0.49999);
    
    % Uffx(((z_global-ht).^2+x_global.^2).^0.5<=Rt)=0;
    % Uffz(((z_global-ht).^2+x_global.^2).^0.5<=Rt)=0;
% %     Uffx(((z_global-ht).^2+x_global.^2).^0.5<=Rt)=Uffx_crown;
% %     Uffz(((z_global-ht).^2+x_global.^2).^0.5<=Rt)=Uffz_crown;
    Uffx(((z_global-ht).^2+x_global.^2).^0.5<Rt)=Uffx_crown;
    Uffz(((z_global-ht).^2+x_global.^2).^0.5<Rt)=Uffz_crown;
    
elseif switch_soil==9000
    [ uxdeep1,uzdeep1 ]=u_verr_deep( z_global,x_global,ht,epsilon,rho_deep,nis,switch_shape);
    Uffx=uxdeep1;
    Uffz =uzdeep1;
    Uffx=coeffux*Uffx;
    Uffz=coeffuz*Uffz;
    Uffy=zeros(nnodes,npali);

    
elseif switch_soil==3
    Uffx=Sx;
    Uffz =Sz;
    Uffx=coeffux*Uffx;
    Uffz=coeffuz*Uffz;
    Uffy=zeros(nnodes,npali);

    
end


Sx=Uffx;
Sy=Uffy;
Sz=Uffz;



%==========================================================================
%CARATTERISTICHE DELL'INPUT (SPOSTAMENTO DEL VINCOLO:"CEDIMENTO VINCOLARE")
%==========================================================================
%Quanto segue definisce l'input degli spostamenti nodali free-field. Il
%programma definisce la matrice Dfft che ha tante colonne per quanti sono
%gli istanti di analisi (comprese code di zeri) e tante righe quanti sono i
%gradi di libertà del singolo palo. I gradi di libertà sono ordinati a
%partire dal nodo in testa al palo secondo 
%  - spostamento orizzontale lungo x
%  - spostamento orizzontale lungo y
%  - spostamento orizzontale lungo z
%  - rotazione intorno a x
%  - rotazione intorno a y
%  - rotazione intorno a z

  
%Uff VECTOR AT NODES OF FINITE ELEMENTS
Dfft_gruppo=zeros(npali*nnodes*6,1);
for i=1:nnodes;
   for ip=1:npali;
    Dfft_gruppo((ip-1)*(nnodes*6)+(i-1)*6+1,1)=Uffx(i,ip);
    Dfft_gruppo((ip-1)*(nnodes*6)+(i-1)*6+2,1)=Uffy(i,ip);
    Dfft_gruppo((ip-1)*(nnodes*6)+(i-1)*6+3,1)=Uffz(i,ip);
   end
end




%% TUNNEL-PILE GROUP IN LAYERED SOIL
% INTERACTION LIMITED WITHIN LAYERS

% STIFFNESS MATRIX OF PILE GROUP
KKpalo=zeros(npali*(length(dz)*6+6),npali*(length(dz)*6+6));
Kpalo=zeros(length(dz)*6+6,length(dz)*6+6,npali,npali);

for i=1:npali
    
    for ii=1:length(dz);  
    KBern3Delt =KBern3D(Ep,d,dz(ii),alfa(i),beta(i));

    %Original code no clash
    % Kpalo((ii-1)*6+1:(ii-1)*6+12,(ii-1)*6+1:(ii-1)*6+12,i)=Kpalo((ii-1)*6+1:(ii-1)*6+12,(ii-1)*6+1:(ii-1)*6+12,i)+KBern3Delt;
    
    

    %Code for clash
    clash_coeff=10000;
    KBern3Delt_clash=KBern3D(Ep/clash_coeff,d,dz(ii),alfa(i),beta(i));
    if x_global(1,i)==0 &&  (sum(dz(1:ii))>ht-Rt)
    Kpalo((ii-1)*6+1:(ii-1)*6+12,(ii-1)*6+1:(ii-1)*6+12,i)=Kpalo((ii-1)*6+1:(ii-1)*6+12,(ii-1)*6+1:(ii-1)*6+12,i)+KBern3Delt_clash;
    else 
    Kpalo((ii-1)*6+1:(ii-1)*6+12,(ii-1)*6+1:(ii-1)*6+12,i)=Kpalo((ii-1)*6+1:(ii-1)*6+12,(ii-1)*6+1:(ii-1)*6+12,i)+KBern3Delt;
    end

    
    
    end
    
    KKpalo((i-1)*((length(dz)*6+6))+1:(i)*((length(dz)*6+6)),(i-1)*((length(dz)*6+6))+1:(i)*((length(dz)*6+6)))=...
    KKpalo((i-1)*((length(dz)*6+6))+1:(i)*((length(dz)*6+6)),(i-1)*((length(dz)*6+6))+1:(i)*((length(dz)*6+6)))+...
    Kpalo(:,:,i);

end

%% Pile stiffness at the pile heads to be reduced 


% STIFFNESS MATRIX OF PILE GROUP
KKpalo_negative=zeros(npali*(length(dz)*6+6),npali*(length(dz)*6+6));
Kpalo_negative=zeros(length(dz)*6+6,length(dz)*6+6,npali,npali);

for i=1:npali
    
    for ii=1:1%length(dz) 
    KBern3Delt =KBern3D(Ep,d,dz(ii),alfa(i),beta(i));

    KBern3Delt_axial=0*KBern3Delt;
    KBern3Delt_axial(:,3)=KBern3Delt(:,3);
    KBern3Delt_axial(:,9)=KBern3Delt(:,9);

    
    Kpalo_negative((ii-1)*6+1:(ii-1)*6+12,(ii-1)*6+1:(ii-1)*6+12,i)=Kpalo_negative((ii-1)*6+1:(ii-1)*6+12,(ii-1)*6+1:(ii-1)*6+12,i)+KBern3Delt_axial;

    
    end
    
    KKpalo_negative((i-1)*((length(dz)*6+6))+1:(i)*((length(dz)*6+6)),(i-1)*((length(dz)*6+6))+1:(i)*((length(dz)*6+6)))=...
    KKpalo_negative((i-1)*((length(dz)*6+6))+1:(i)*((length(dz)*6+6)),(i-1)*((length(dz)*6+6))+1:(i)*((length(dz)*6+6)))+...
    Kpalo_negative(:,:,i);

end

%% Stiffness matrix to simulate very stiff cap
KKcap_fixed=zeros(npali*(length(dz)*6+6),npali*(length(dz)*6+6));
Kcap_fixed=zeros(npali*(6),npali*(6));
% 

% %segmetal beam connecting pile heads
% % %Stiffness matrix only considering the pile head nodes
% % for ii=1:npali-1 
% %     Xi=[X(ii)   Y(ii) 0];
% %     Xf=[X(ii+1) Y(ii+1) 0 ];
% %     [Kelem_cap] = KBern3D_CAP(Ep*100000,d,Xi,Xf);
% %   
% %     Kcap_fixed((ii-1)*6+1:(ii-1)*6+12,(ii-1)*6+1:(ii-1)*6+12)=...
% %     Kcap_fixed((ii-1)*6+1:(ii-1)*6+12,(ii-1)*6+1:(ii-1)*6+12)+Kelem_cap;
% % 
% % end
% % 
% % for ii=1:npali
% %     for jj=1:npali
% %     KKcap_fixed((ii-1)*nnodes*6+1:(ii-1)*nnodes*6+6,(jj-1)*nnodes*6+1:(jj-1)*nnodes*6+6)=...
% %     Kcap_fixed((ii-1)*6+1:(ii-1)*6+6,(jj-1)*6+1:(jj-1)*6+6);
% %     end
% % end

%all pile heads connected to others
%Stiffness matrix only considering the pile head nodes
for ii=1:npali 
    
for jj=1:npali
    
    if ii==jj
    else
        
    Xi=[X(ii)   Y(ii) 0];
    Xf=[X(jj)   Y(jj) 0 ];
    [Kelem_cap] = KBern3D_CAP(Ep*1,d*3,Xi,Xf);%d*3
  
    KKcap_fixed((ii-1)*nnodes*6+1:(ii-1)*nnodes*6+6,(ii-1)*nnodes*6+1:(ii-1)*nnodes*6+6)=...
    KKcap_fixed((ii-1)*nnodes*6+1:(ii-1)*nnodes*6+6,(ii-1)*nnodes*6+1:(ii-1)*nnodes*6+6)+...
    Kelem_cap(1:6,1:6);
    
    KKcap_fixed((ii-1)*nnodes*6+1:(ii-1)*nnodes*6+6,(jj-1)*nnodes*6+1:(jj-1)*nnodes*6+6)=...
    KKcap_fixed((ii-1)*nnodes*6+1:(ii-1)*nnodes*6+6,(jj-1)*nnodes*6+1:(jj-1)*nnodes*6+6)+...
    Kelem_cap(1:6,7:12);

    KKcap_fixed((jj-1)*nnodes*6+1:(jj-1)*nnodes*6+6,(jj-1)*nnodes*6+1:(jj-1)*nnodes*6+6)=...
    KKcap_fixed((jj-1)*nnodes*6+1:(jj-1)*nnodes*6+6,(jj-1)*nnodes*6+1:(jj-1)*nnodes*6+6)+...
    Kelem_cap(7:12,7:12);
    
    KKcap_fixed((jj-1)*nnodes*6+1:(jj-1)*nnodes*6+6,(ii-1)*nnodes*6+1:(ii-1)*nnodes*6+6)=...
    KKcap_fixed((jj-1)*nnodes*6+1:(jj-1)*nnodes*6+6,(ii-1)*nnodes*6+1:(ii-1)*nnodes*6+6)+...
    Kelem_cap(7:12,1:6);
    end
end
end

if switch_cap==2
KKframe_fixed =KKframe_fixed  +KKcap_fixed;
KKframe_hinged=KKframe_hinged +KKcap_fixed;
end

%% test
% % force=[11309.9848803355]/2;
% % % force/Ep/(pi*d^2/4)*L;
% % 
% % F(3)=force;
% % F(39)=force;
% % F(72)=0;
% % KK=KKpalo*0;
% % val=10^15;
% % KK(36,36)=val;
% % KK(35,35)=val;
% % KK(34,34)=val;
% % KK(33,33)=val;
% % KK(32,32)=val;
% % KK(31,31)=val;
% % 
% % KK(72,72)=val;
% % KK(71,71)=val;
% % KK(70,70)=val;
% % KK(69,69)=val;
% % KK(68,68)=val;
% % KK(67,67)=val;
% % 
% % beta=0.95;
% % KKnew=KKpalo-beta*KKpalo_negative+KK;
% % U=inv(KKnew)*F';
% % (U(3)-U(33))
% % (U(3)-U(33))/((force/Ep/(pi*d^2/4)*(L-2))+(force/((1-beta)*Ep)/(pi*d^2/4)*(2)))
% % 
% % ((force/((1-0.9)*Ep)/(pi*d^2/4)*(2)))


%%

%
% 
% =========================================================================
% CALCOLO DELLA MATRICE DI VINCOLO A
% =========================================================================

%d=AdE
%d=[s1,...,sp,...,sn]t   dove sp spostamenti nodali del pali p-esimo
%dE[sE1,...,sEp,...,sEn]t   dove sEp spostamenti nodali del pali p-esimo a
%meno degli spostamenti a quota z=0

A_cap_i=zeros((6*length(dz)+6)*npali,6+(6*length(dz))*npali);
A_cap_h=zeros((6*length(dz)+6)*npali,6+(6*length(dz)+3)*npali);

for i=1:npali;
    A_cap_i((i-1)*(6*length(dz)+6)+1:(i-1)*(6*length(dz)+6)+6,1:6)=Ai(Xm,Ym,Zm,X(i),Y(i),z(1));   
    A_cap_i(7+(i-1)*(6*length(dz)+6):(i)*(6*length(dz)+6),7+(i-1)*(6*length(dz)):6+(i)*(6*length(dz)))=eye(6*length(dz),6*length(dz));
end

for i=1:npali;
    A_cap_h((i-1)*(6*length(dz)+6)+1:(i-1)*(6*length(dz)+6)+3,1:6)=Ah(Xm,Ym,Zm,X(i),Y(i),z(1));   
    A_cap_h(4+(i-1)*(6*length(dz)+6):(i)*(6*length(dz)+6),7+(i-1)*(6*length(dz)+3):6+(i)*(6*length(dz)+3))=eye(6*length(dz)+3,6*length(dz)+3);
end









%% Kitiyodom METHOD
% SPRINGS INTERACTION ACCORDING TO MINDLIN'S SOLUTION

% KKsoil_jap_integration=zeros(npali*(length(dz)*6+6),npali*(length(dz)*6+6));
KKsoil=zeros(npali*(nelementi+1)*6);
KKsoil_base=zeros(npali*(nelementi+1)*6);
% 
%Winkler model of the soil
rm=2.5*L*(1-nis);
[ kh,kv ] = soilsprings_static(rm,Es(1),Ep(1),d,nis);
kb=d*Eb/(1-nib^2);



FLEX=zeros(npali*(nelementi+1)*3);

if  exist('switch_type')
    if switch_type==0
        nint=8*4;
        %FLEXIBILITY MATRIX for the continuum
        for jg=1:npali*(nelementi+1)

            for ii_int=1:nint



        %     Zj_nod(1:1:npali*(nelementi+1),1)=Zi_nod(jg);
        %     Xj_nod(1:1:npali*(nelementi+1),1)=Xi_nod(jg)+d/2*cos(2*pi*(ii_int-1)/nint);
        %     Yj_nod(1:1:npali*(nelementi+1),1)=Yi_nod(jg)+d/2*sin(2*pi*(ii_int-1)/nint);
        % 


            if ii_int<=nint/2

            Zj_nod(1:1:npali*(nelementi+1),1)=Zi_nod(jg);
            Xj_nod(1:1:npali*(nelementi+1),1)=Xi_nod(jg)+d/2*cos(2*pi*(ii_int-1)/nint*2);
            Yj_nod(1:1:npali*(nelementi+1),1)=Yi_nod(jg)+d/2*sin(2*pi*(ii_int-1)/nint*2);

            else
            Zj_nod(1:1:npali*(nelementi+1),1)=Zi_nod(jg);
            Xj_nod(1:1:npali*(nelementi+1),1)=Xi_nod(jg)+d/4*cos(2*pi*((ii_int-nint/2)-1)/nint*2);
            Yj_nod(1:1:npali*(nelementi+1),1)=Yi_nod(jg)+d/4*sin(2*pi*((ii_int-nint/2)-1)/nint*2);

            end    

            dZj_nod(:,jg)=Zi_nod(:)-Zj_nod(:);
            dXj_nod(:,jg)=Xi_nod(:)-Xj_nod(:);
            dYj_nod(:,jg)=Yi_nod(:)-Yj_nod(:); 







            [ uxij_1,uyij_1,uzij_1 ] = int_factor_mindlin_j_1_vect_CONT( Xi_nod,Yi_nod,Zi_nod,Xj_nod,Yj_nod,Zj_nod,Gs_nod,ni_nod);
            [ uxij_2,uyij_2,uzij_2 ] = int_factor_mindlin_j_2_vect_CONT( Xi_nod,Yi_nod,Zi_nod,Xj_nod,Yj_nod,Zj_nod,Gs_nod,ni_nod);
            [ uxij_3,uyij_3,uzij_3 ] = int_factor_mindlin_j_3_vect_CONT( Xi_nod,Yi_nod,Zi_nod,Xj_nod,Yj_nod,Zj_nod,Gs_nod,ni_nod);





            FLEX(1:3:npali*(nelementi+1)*3,(jg-1)*3+1)=uxij_1/nint+FLEX(1:3:npali*(nelementi+1)*3,(jg-1)*3+1);
            FLEX(2:3:npali*(nelementi+1)*3,(jg-1)*3+1)=uyij_1/nint+FLEX(2:3:npali*(nelementi+1)*3,(jg-1)*3+1);
            FLEX(3:3:npali*(nelementi+1)*3,(jg-1)*3+1)=uzij_1/nint+FLEX(3:3:npali*(nelementi+1)*3,(jg-1)*3+1);

            FLEX(1:3:npali*(nelementi+1)*3,(jg-1)*3+2)=uxij_2/nint+FLEX(1:3:npali*(nelementi+1)*3,(jg-1)*3+2);
            FLEX(2:3:npali*(nelementi+1)*3,(jg-1)*3+2)=uyij_2/nint+FLEX(2:3:npali*(nelementi+1)*3,(jg-1)*3+2);
            FLEX(3:3:npali*(nelementi+1)*3,(jg-1)*3+2)=uzij_2/nint+FLEX(3:3:npali*(nelementi+1)*3,(jg-1)*3+2);

            FLEX(1:3:npali*(nelementi+1)*3,(jg-1)*3+3)=uxij_3/nint+FLEX(1:3:npali*(nelementi+1)*3,(jg-1)*3+3);
            FLEX(2:3:npali*(nelementi+1)*3,(jg-1)*3+3)=uyij_3/nint+FLEX(2:3:npali*(nelementi+1)*3,(jg-1)*3+3);
            FLEX(3:3:npali*(nelementi+1)*3,(jg-1)*3+3)=uzij_3/nint+FLEX(3:3:npali*(nelementi+1)*3,(jg-1)*3+3);



            end




        end
        
    elseif switch_type==1
    run('Main_2');
    end
end

if  exist('switch_chenflex')

   if switch_chenflex==1   
    load FLEX_CHEN_Huang_Singlepile-singlelayer.mat
   elseif switch_chenflex==2
    load FLEX_CHEN_Huang_Singlepile-twolayer.mat   
    
   elseif switch_chenflex==11
    load Case1.mat 
    elseif switch_chenflex==15
    load Case5.mat   
    elseif switch_chenflex==16
    load Case6.mat    

    
    
    elseif switch_chenflex==221
    load case221.mat ; FLEX=FLEX221; 
    elseif switch_chenflex==222
    load case222.mat ; FLEX=FLEX222;
    elseif switch_chenflex==223
    load case223.mat ; FLEX=FLEX223;
    elseif switch_chenflex==224
    load case224.mat ; FLEX=FLEX224;
    elseif switch_chenflex==241
    load case241.mat ; FLEX=FLEX241;
    elseif switch_chenflex==242
    load case242.mat ; FLEX=FLEX242;
    elseif switch_chenflex==243
    load case243.mat ; FLEX=FLEX243;
    elseif switch_chenflex==244
    load case244.mat ; FLEX=FLEX244;
    
    
    elseif switch_chenflex==500
    load FLEX_BrianPG_3x3_s50d.mat ; FLEX=FLEX; 
    elseif switch_chenflex==499
    load FLEX_BrianPG_3x3_s25d.mat ; FLEX=FLEX; 
    elseif switch_chenflex==501
    load FLEX_BrianSP_3.mat ; FLEX=FLEX; 
    
   elseif switch_chenflex==1000
    load Case0_PG.mat ; FLEX=FLEX;   
    elseif switch_chenflex==1001
    load Case1_PG.mat ; FLEX=FLEX;  
    elseif switch_chenflex==1002
    load Case2_PG.mat ; FLEX=FLEX;  
    
    elseif switch_chenflex==20000
    load matlab_flex_br.mat;  
    
    elseif switch_chenflex==600
    load Flex_Ong.mat;  
    
        
   end
    
end

% %FLEXIBILITY MATRIX WITHOUT INTERACTION BETWEEN NODES IN THE SAME PILE
% %Interection effects between the nodes within the same pile is neglected:
% %components outside diagonal are set to zero if nodes i and j with same
% %pile
% 
% % for in=1:(npali*(nelementi+1)*3)
% %     for jn=1:(npali*(nelementi+1)*3)
% %         i=ceil(in/(nelementi+1)/3);  %ie
% %         j=ceil(jn/(nelementi+1)/3);  %je
% %         if i==j &&  in~=jn
% %         FLEX(in,jn)=0;
% %         end
% %     end
% % end



%STIFFNESS MATRIX OF THE SOIL
%Only translational dofs
KKsoil_rid=inv(FLEX);


%Zeros are added within soil stiffness matrix to account for rotational
%dofs
for i=1:npali*(nelementi+1)
    for j=1:npali*(nelementi+1)
       KKsoil(1+(i-1)*6:3+(i-1)*6,1+(j-1)*6:3+(j-1)*6)=...
        KKsoil_rid(1+(i-1)*3:3+(i-1)*3,1+(j-1)*3:3+(j-1)*3);
	end
end



%TORSIONAL SPRINGS AT PILE TIP TO REMOVE UNCOSTRAINED TORSIONAL DEGREE OF FREEDOM 
KKtors=zeros(6*nnodes*npali);
for i=1:npali
	for j=1:npali
        if i==j 
            KKtors((i-1)*((nnodes*6))+6+6*(nnodes-1),(i-1)*((nnodes*6))+6+6*(nnodes-1))=max(max(KKsoil));%1
        end
	end
end


%GLOBAL STIFFNESS MATRIX    
    KKpg=KKpalo+KKsoil+KKtors;  %KKsoil_base+
    
%GLOBAL STIFFNESS MATRIX FRAME WITH FIXED PILE HEADS  
    KKpg_fr_i=KKframe_fixed+KKpalo+KKsoil+KKtors;  %+KKsoil_base

%GLOBAL STIFFNESS MATRIX FRAME WITH HINGED PILE HEADS  
    KKpg_fr_h=KKframe_hinged+KKpalo+KKsoil+KKtors; %+KKsoil_base
    
%GLOBAL STIFFNESS MATRIX CONSTRAINED BY RIGID CAP WITH FIXED HEADS
    KKpg_cap_i=A_cap_i'*KKpg*A_cap_i;
    
%GLOBAL STIFFNESS MATRIX CONSTRAINED BY RIGID CAP WITH HINGED HEADS
    KKpg_cap_h=A_cap_h'*(KKpg)*A_cap_h; 

    if length(X)==1
        KKpg_cap_h(4,4)= KKpg_cap_h(4,4)/10^10;
    end
    
%EXTERNAL FORCES
Fpg=(KKsoil)*Dfft_gruppo; %vettore forze esterne %+KKsoil_base

%EXTERNAL FORCES AT THE PILE HEAD
% Code to inser unit force at pile head
% Fpg=Fpg*0;
% P_pilehead=[450*10^3];
if exist('P_pilehead')
    for i=1:npali
        Fpg(3+6*nnodes*(i-1),1)=Fpg(3+6*nnodes*(i-1),1)+P_pilehead;
    end
end
% M_pilehead=P_pilehead*9*d*5;
if exist('M_pilehead')
        Fpg(5,1)=Fpg(5,1)+M_pilehead;
end

%EXTERNAL FORCES IN CASE OF RIGID CAP WITH FIXED HEADS
Fpg_cap_i=A_cap_i'*Fpg;             %vettore forze esterne nel dominio delle frequenze dopo la condensazione

%EXTERNAL FORCES IN CASE OF RIGID CAP WITH HINGED HEADS
Fpg_cap_h=A_cap_h'*Fpg;             %vettore forze esterne nel dominio delle frequenze dopo la condensazione

%SOLUTION
U_cap_i=(KKpg_cap_i)\Fpg_cap_i;  %Foundation input motion=campo di spostamenti soluzione condensato (nell'articolo d-tilde)
U_node_cap_i=A_cap_i*U_cap_i(:,1);                      %campo di spostamenti soluzione di tutti i nodi (compresi i nodi slave) (nell'articolo d)
U_cap_i_master=[U_cap_i(3,1),U_cap_i(1,1),U_cap_i(5,1)];

%SOLUTION
U_cap_h=(KKpg_cap_h)\Fpg_cap_h;  %Foundation input motion=campo di spostamenti soluzione condensato (nell'articolo d-tilde)
U_node_cap_h=A_cap_h*U_cap_h(:,1);                      %campo di spostamenti soluzione di tutti i nodi (compresi i nodi slave) (nell'articolo d)
U_cap_h_master=[U_cap_h(3,1),U_cap_h(1,1),U_cap_h(5,1)];


%SOLUTION     
U_node_fr_i=(KKpg_fr_i)\Fpg;
U_node_fr_h=(KKpg_fr_h)\Fpg; 
U_node_free=(KKpg)\Fpg; 

   



for i=1:nnodes;
	for ip=1:npali;
	Fxpg(i,ip)=Fpg((ip-1)*(nnodes*6)+((i-1)*6)+1);
    Fypg(i,ip)=Fpg((ip-1)*(nnodes*6)+((i-1)*6)+2);
    Fzpg(i,ip)=Fpg((ip-1)*(nnodes*6)+((i-1)*6)+3);
    
    Uxpg_cap_i(i,ip)=U_node_cap_i((ip-1)*(nnodes*6)+((i-1)*6)+1);
    Uypg_cap_i(i,ip)=U_node_cap_i((ip-1)*(nnodes*6)+((i-1)*6)+2);
    Uzpg_cap_i(i,ip)=U_node_cap_i((ip-1)*(nnodes*6)+((i-1)*6)+3);
    Ur2pg_cap_i(i,ip)=U_node_cap_i((ip-1)*(nnodes*6)+((i-1)*6)+5);

    Uxpg_cap_h(i,ip)=U_node_cap_h((ip-1)*(nnodes*6)+((i-1)*6)+1);
    Uypg_cap_h(i,ip)=U_node_cap_h((ip-1)*(nnodes*6)+((i-1)*6)+2);
    Uzpg_cap_h(i,ip)=U_node_cap_h((ip-1)*(nnodes*6)+((i-1)*6)+3);
    Ur2pg_cap_h(i,ip)=U_node_cap_h((ip-1)*(nnodes*6)+((i-1)*6)+5);
    
    Uxpg_fr_i(i,ip)=U_node_fr_i((ip-1)*(nnodes*6)+((i-1)*6)+1);
    Uypg_fr_i(i,ip)=U_node_fr_i((ip-1)*(nnodes*6)+((i-1)*6)+2);
    Uzpg_fr_i(i,ip)=U_node_fr_i((ip-1)*(nnodes*6)+((i-1)*6)+3);
    Ur2pg_fr_i(i,ip)=U_node_fr_i((ip-1)*(nnodes*6)+((i-1)*6)+5);
    
    Uxpg_fr_h(i,ip)=U_node_fr_h((ip-1)*(nnodes*6)+((i-1)*6)+1);
    Uypg_fr_h(i,ip)=U_node_fr_h((ip-1)*(nnodes*6)+((i-1)*6)+2);
    Uzpg_fr_h(i,ip)=U_node_fr_h((ip-1)*(nnodes*6)+((i-1)*6)+3);   
    Ur2pg_fr_h(i,ip)=U_node_fr_h((ip-1)*(nnodes*6)+((i-1)*6)+5);
    
    Uxpg_free(i,ip)=U_node_free((ip-1)*(nnodes*6)+((i-1)*6)+1);
    Uypg_free(i,ip)=U_node_free((ip-1)*(nnodes*6)+((i-1)*6)+2);
    Uzpg_free(i,ip)=U_node_free((ip-1)*(nnodes*6)+((i-1)*6)+3);
    Ur2pg_free(i,ip)=U_node_free((ip-1)*(nnodes*6)+((i-1)*6)+5);
    
    end
end

%INNER FORCES
Ip=pi*(d/2)^4/4;
for i_p=1:npali
    for ii=1:nnodes
        if ii==1
            Mpg_cap_i(ii,i_p)=-Ep*Ip*(Ur2pg_cap_i(ii,i_p)-Ur2pg_cap_i(ii+1,i_p))/h_el;
            Mpg_cap_h(ii,i_p)=-Ep*Ip*(Ur2pg_cap_h(ii,i_p)-Ur2pg_cap_h(ii+1,i_p))/h_el;
            Mpg_free(ii,i_p)=Ep*Ip*(Uxpg_free(ii,i_p)-2*Uxpg_free(ii+1,i_p)+Uxpg_free(ii+2,i_p))/h_el^2;
            Mpg_fr_h(ii,i_p)=KKframe_hinged(5+6*(nnodes)*(i_p-1),:)*U_node_fr_h;
            Mpg_fr_i(ii,i_p)=KKframe_fixed(5+6*(nnodes)*(i_p-1),:)*U_node_fr_i;
        elseif ii==nnodes
            Mpg_cap_i(ii,i_p)=0;  
            Mpg_cap_h(ii,i_p)=0;
            Mpg_free(ii,i_p)=0;
            Mpg_fr_h(ii,i_p)=0;
            Mpg_fr_i(ii,i_p)=0;  
        else
            Mpg_cap_i(ii,i_p)=Ep*Ip*(Uxpg_cap_i(ii-1,i_p)-2*Uxpg_cap_i(ii,i_p)+Uxpg_cap_i(ii+1,i_p))/h_el^2;
            Mpg_cap_h(ii,i_p)=Ep*Ip*(Uxpg_cap_h(ii-1,i_p)-2*Uxpg_cap_h(ii,i_p)+Uxpg_cap_h(ii+1,i_p))/h_el^2;
            Mpg_free(ii,i_p)=Ep*Ip*(Uxpg_free(ii-1,i_p)-2*Uxpg_free(ii,i_p)+Uxpg_free(ii+1,i_p))/h_el^2;
            Mpg_fr_h(ii,i_p)=Ep*Ip*(Uxpg_fr_h(ii-1,i_p)-2*Uxpg_fr_h(ii,i_p)+Uxpg_fr_h(ii+1,i_p))/h_el^2;
            Mpg_fr_i(ii,i_p)=Ep*Ip*(Uxpg_fr_i(ii-1,i_p)-2*Uxpg_fr_i(ii,i_p)+Uxpg_fr_i(ii+1,i_p))/h_el^2;
        end
    end
end

for i_p=1:npali
    for ii=1:nnodes
        if ii==1
            Npg_cap_i(ii,i_p)=Ep(1)*Ap*(Uzpg_cap_i(ii+1,i_p)-Uzpg_cap_i(ii,i_p))/h_el;
            Npg_cap_h(ii,i_p)=Ep(1)*Ap*(Uzpg_cap_h(ii+1,i_p)-Uzpg_cap_h(ii,i_p))/h_el;
            Npg_free(ii,i_p)=Ep(1)*Ap*(Uzpg_free(ii+1,i_p)-Uzpg_free(ii,i_p))/h_el;
            Npg_fr_h(ii,i_p)=KKframe_hinged(3+6*(nnodes)*(i_p-1),:)*U_node_fr_h;
            Npg_fr_i(ii,i_p)=KKframe_fixed(3+6*(nnodes)*(i_p-1),:)*U_node_fr_i;            %F0 ext loads RR_hinhead*UU_head frame loads
        elseif ii==nnodes          
            Npg_cap_i(ii,i_p)=Ep(1)*Ap*(Uzpg_cap_i(ii,i_p)-Uzpg_cap_i(ii-1,i_p))/h_el; 
            Npg_cap_h(ii,i_p)=Ep(1)*Ap*(Uzpg_cap_h(ii,i_p)-Uzpg_cap_h(ii-1,i_p))/h_el;
            Npg_free(ii,i_p)=Ep(1)*Ap*(Uzpg_free(ii,i_p)-Uzpg_free(ii-1,i_p))/h_el; 
            Npg_fr_h(ii,i_p)=Ep*Ap*(Uzpg_fr_h(ii,i_p)-Uzpg_fr_h(ii-1,i_p))/h_el;
            Npg_fr_i(ii,i_p)=Ep*Ap*(Uzpg_fr_i(ii,i_p)-Uzpg_fr_i(ii-1,i_p))/h_el;    %note the Fztpg_spr_int
        else
            Npg_cap_i(ii,i_p)=Ep(1)*Ap*(Uzpg_cap_i(ii+1,i_p)-Uzpg_cap_i(ii-1,i_p))/2/h_el;
            Npg_cap_h(ii,i_p)=Ep(1)*Ap*(Uzpg_cap_h(ii+1,i_p)-Uzpg_cap_h(ii-1,i_p))/2/h_el;
            Npg_free(ii,i_p)=Ep(1)*Ap*(Uzpg_free(ii+1,i_p)-Uzpg_free(ii-1,i_p))/2/h_el;
            Npg_fr_h(ii,i_p)=Ep*Ap*(Uzpg_fr_h(ii+1,i_p)-Uzpg_fr_h(ii-1,i_p))/2/h_el;
            Npg_fr_i(ii,i_p)=Ep*Ap*(Uzpg_fr_i(ii+1,i_p)-Uzpg_fr_i(ii-1,i_p))/2/h_el;
        end
    end
end




%% SIMPLIFIED METHODS FOR TUNNEL SINGLE PILE SETTLEMENTS
%PILE SETTLEMENT = 2/3 GF

Uzpg_twothird=zeros(1,npali);
for i_p=1:npali
    i=fix(2/3*(nnodes));
    Uzpg_twothird(1,i_p)=Sz(i,i_p);  

end

%PILE SETTLEMENT = AVG GF
Uzpg_avg=zeros(1,npali);
for i_p=1:npali
    for i=1:nnodes
        Uzpg_avg(1,i_p)=Sz(i,i_p)/(nnodes)+Uzpg_avg(1,i_p);  
    end
end


%PILE SETTLEMENT FOR RIGID PILE
Uzpg_rig_base=zeros(1,npali);
Uzpg_rig_base_num=zeros(1,npali);
Uzpg_rig_base_den=zeros(1,npali);
for i_p=1:npali
    for i=1:nnodes
        if i==nnodes
         Uzpg_rig_base_num(1,i_p)=Sz(i,i_p)*(h_el/2*kv+kb)+Uzpg_rig_base_num(1,i_p); 
        elseif i==1
         Uzpg_rig_base_num(1,i_p)=Sz(i,i_p)*(h_el/2*kv+0)+Uzpg_rig_base_num(1,i_p); 
        else
         Uzpg_rig_base_num(1,i_p)=Sz(i,i_p)*(h_el*kv+0)+Uzpg_rig_base_num(1,i_p); 
        end
    end
    Uzpg_rig_base_den(1,i_p)=(nnodes-1)*h_el*kv+kb;
    Uzpg_rig_base=Uzpg_rig_base_num./Uzpg_rig_base_den;
end

%%
% =========================================================================
% INPUT
% =========================================================================

nfoot=npali;
nnodes_foot=nnodes;


% PARAMETERS
    if switch_nep==1
        %vertical load
%         num_incr_TUNNEL=200;
%         num_incr_P=4*50;
%         delta_err   =1.0e-2; 
%         eps_err     =1.0e-2; 
%         beta=1.0;




%         tunnelling
        num_incr_TUNNEL=500;        %500 
        num_incr_P=4*60;            %4*60
        delta_err   =1.0e-2; 
        eps_err     =1.0e-2; 
        beta=1.0;
        
        if coeffux==0
            if coeffuz==0
            num_incr_TUNNEL=100;        %500 
            num_incr_P=100;            %4*60
            beta=3;
            end
        end
        
    else
        num_incr_TUNNEL=50;
        num_incr_P=20;%20
        delta_err   =1.0e-3; %set at least at 1.0e-3
        eps_err     =1.0e-3; %set at least at 1.0e-3
        beta=1.0;
        
        
        
        if coeffux==0
            if coeffuz==0
            num_incr_TUNNEL=100;        %500 
            num_incr_P=100;            %4*60
            beta=3;
            end
        end
        
        
    end



% flim2=0;
% flim1=Inf;
% mu=1000*tan(2*pi/360*phi_int);    %friction coefficient

Vltp_vect=0:Vltp/num_incr_TUNNEL:Vltp;
epsilon_vect=Vltp_vect/100/2;

%% LIMIT FORCES

if switch_limf==1 %Alpha SOLUTION
    


%Vertical direction
cu=cu0+cuslope.*z';
flim2_vector       =-d*pi*h_el.*alpha_cu.*cu;%Compression         %maximum tensile force of plastic elements (<=0)
flim1_vector       =+d*pi*h_el.*alpha_cu.*cu;%Tensile Inf;%It needs to be set as a limit prssure rather than a force considering area corresponding to each node%2.5*10^5/(0.5*1.2);%       %maximum compressive force of plastic elements   (>=0) %     flim1       = c *Nc *sc +  0.5 gamma * B *Ng *sg;
flim2_vector(1)=flim2_vector(1)/2;
flim1_vector(1)=flim1_vector(1)/2;
flim2_vector(end)=flim2_vector(end)/2;
flim1_vector(end)=qb_lim_alpha*pi*d^2/4+flim1_vector(end)/2;%/2
    flim2_vector_singlepile=flim2_vector;
    flim1_vector_singlepile=flim1_vector;
    flim2_vector=repmat(flim2_vector,npali,1);
    flim1_vector=repmat(flim1_vector,npali,1);
        flim2_vector_vlt=repmat(flim2_vector,1,length(Vltp_vect));
        flim1_vector_vlt=repmat(flim1_vector,1,length(Vltp_vect));

%Horizontal direction
xlim2_vectora       =-d*h_el.*9.*cu;
xlim1_vectora       =+d*h_el.*9.*cu;
xlim2_vectorb       =-d*h_el.*2*(1+z'/d).*cu;
xlim1_vectorb       =+d*h_el.*2*(1+z'/d).*cu;
% xlim2_vectorb       =-d*h_el.*(2+2.333.*z'/d).*cu;
% xlim1_vectorb       =+d*h_el.*(2+2.33.*z'/d).*cu;

if exist('switch_nohorizontalalphamethod')==1
%Horizontal direction
xlim2_vectora       =-d*h_el.*9.*cu/10^3;
xlim1_vectora       =+d*h_el.*9.*cu/10^3;
xlim2_vectorb       =-d*h_el.*2*(1+z'/d).*cu/10^3;
xlim1_vectorb       =+d*h_el.*2*(1+z'/d).*cu/10^3;
% xlim2_vectorb       =-d*h_el.*(2+2.333.*z'/d).*cu;
% xlim1_vectorb       =+d*h_el.*(2+2.33.*z'/d).*cu; 
    
end

xlim2_vector=max(xlim2_vectora,xlim2_vectorb);
xlim1_vector=min(xlim1_vectora,xlim1_vectorb);

xlim2_vector(1)=xlim2_vector(1)/2;
xlim1_vector(1)=xlim1_vector(1)/2;
xlim2_vector(end)=xlim2_vector(end)/2;
xlim1_vector(end)=xlim1_vector(end)/2;
    xlim2_vector_singlepile=xlim2_vector;
    xlim1_vector_singlepile=xlim1_vector;
    xlim2_vector=repmat(xlim2_vector,npali,1);
    xlim1_vector=repmat(xlim1_vector,npali,1);
        xlim2_vector_vlt=repmat(xlim2_vector,1,length(Vltp_vect));
        xlim1_vector_vlt=repmat(xlim1_vector,1,length(Vltp_vect));



elseif switch_limf==2 %Beta SOLUTION

sigmavop=k0_coeff.*gammap.*z';
flim2_vector       =-d*pi*h_el.*(beta_flim.*sigmavop+c_beta);%Tensile         %maximum tensile force of plastic elements (<=0)
flim1_vector       =+d*pi*h_el.*(beta_flim.*sigmavop+c_beta);%Compression Inf;%It needs to be set as a limit prssure rather than a force considering area corresponding to each node%2.5*10^5/(0.5*1.2);%       %maximum compressive force of plastic elements   (>=0) %     flim1       = c *Nc *sc +  0.5 gamma * B *Ng *sg;
flim2_vector(1)=flim2_vector(1)/2;
flim1_vector(1)=flim1_vector(1)/2;
flim2_vector(end)=flim2_vector(end)/2;
flim1_vector(end)=qb_lim_beta*pi*d^2/4+flim1_vector(end)/2;%/2/2
    flim2_vector_singlepile=flim2_vector;
    flim1_vector_singlepile=flim1_vector;
    flim2_vector=repmat(flim2_vector,npali,1);
    flim1_vector=repmat(flim1_vector,npali,1);
        flim2_vector_vlt=repmat(flim2_vector,1,length(Vltp_vect));
        flim1_vector_vlt=repmat(flim1_vector,1,length(Vltp_vect));


%Horizontal direction

xlim2_vector       =-d*h_el.*(3*kp_coeff.*sigmavop);
xlim1_vector       =+d*h_el.*(3*kp_coeff.*sigmavop);
xlim2_vector(1)=xlim2_vector(1)/2;
xlim1_vector(1)=xlim1_vector(1)/2;
xlim2_vector(end)=xlim2_vector(end)/2;
xlim1_vector(end)=xlim1_vector(end)/2;
    xlim2_vector_singlepile=xlim2_vector;
    xlim1_vector_singlepile=xlim1_vector;
    xlim2_vector=repmat(xlim2_vector,npali,1);
    xlim1_vector=repmat(xlim1_vector,npali,1);
        xlim2_vector_vlt=repmat(xlim2_vector,1,length(Vltp_vect));
        xlim1_vector_vlt=repmat(xlim1_vector,1,length(Vltp_vect));

elseif switch_limf==4 %tabular data

flim2_vector       =-d*pi*h_el.*zstresslim_vector_tabular;
flim1_vector       =+d*pi*h_el.*zstresslim_vector_tabular;
flim2_vector(1)=flim2_vector(1)/2;
flim1_vector(1)=flim1_vector(1)/2;
flim2_vector(end)=flim2_vector(end)/2;
flim1_vector(end)=qb_lim_alpha*pi*d^2/4+flim1_vector(end)/2;%/2
    flim2_vector_singlepile=flim2_vector;
    flim1_vector_singlepile=flim1_vector;
    flim2_vector=repmat(flim2_vector,npali,1);
    flim1_vector=repmat(flim1_vector,npali,1);
        flim2_vector_vlt=repmat(flim2_vector,1,length(Vltp_vect));
        flim1_vector_vlt=repmat(flim1_vector,1,length(Vltp_vect));

xlim2_vector       =-d*h_el.*xstresslim_vector_tabular;
xlim1_vector       =+d*h_el.*xstresslim_vector_tabular;
xlim2_vector(1)=xlim2_vector(1)/2;
xlim1_vector(1)=xlim1_vector(1)/2;
xlim2_vector(end)=xlim2_vector(end)/2;
xlim1_vector(end)=xlim1_vector(end)/2;
    xlim2_vector_singlepile=xlim2_vector;
    xlim1_vector_singlepile=xlim1_vector;
    xlim2_vector=repmat(xlim2_vector,npali,1);
    xlim1_vector=repmat(xlim1_vector,npali,1);
        xlim2_vector_vlt=repmat(xlim2_vector,1,length(Vltp_vect));
        xlim1_vector_vlt=repmat(xlim1_vector,1,length(Vltp_vect));
    
elseif switch_limf==3 %CE SOLUTION

[blank, SYSname] = system('hostname');
if strcmp(SYSname(1:end-1),'Andrea-UK')
    addpath('C:\Users\Andrea\Dropbox\UoN\MATLAB\AlecMarshall2018CE\')
addpath('C:\Users\Andrea\Dropbox\UoN\MATLAB\AlecMarshall2018CE\CE\')
addpath('C:\Users\Andrea\Dropbox\UoN\MATLAB\AlecMarshall2018CE\Parameters_CE\')

elseif strcmp(SYSname(1:end-1),'MSI')
    addpath('C:\Users\andre\Dropbox (Personal)\UoN\MATLAB\AlecMarshall2018CE\')
addpath('C:\Users\andre\Dropbox (Personal)\UoN\MATLAB\AlecMarshall2018CE\CE\')
addpath('C:\Users\andre\Dropbox (Personal)\UoN\MATLAB\AlecMarshall2018CE\Parameters_CE\')

else


addpath('C:\Users\af624\Dropbox (Personal)\UoN\MATLAB\AlecMarshall2018CE\')
addpath('C:\Users\af624\Dropbox (Personal)\UoN\MATLAB\AlecMarshall2018CE\CE\')
addpath('C:\Users\af624\Dropbox (Personal)\UoN\MATLAB\AlecMarshall2018CE\Parameters_CE\')

end


vltarray=Vltp_vect(2:end); 
% inputarray={'P71N','P61N','P51N','P11N','P21N','P31N','P41N'};
inp.hel=h_el;
fname='f_lim_dataset.mat';
RversusVlt_w(inputarray,fname,vltarray,inp);
load f_lim_dataset.mat


%%
load f_lim_dataset.mat

figure(3)
%set(1, 'Units','centimeters', 'Position',[0 0 18 6])

% a=out_Qbase_0vlt+sum(out_Tshaft_0vlt,1);
% a(1)
% b=+sum(out_Tshaft_0vlt,1);
% b(1)

hold on

for ii=1:size(out_Qbase_0vlt,2)
temp1(:)=out_Qbase_0vlt(:,ii,:);
plot(vltarray_0vlt,temp1,'r','LineWidth',1.5)
temp2(:)=sum(out_Tshaft_0vlt(:,ii,:),1);
plot(vltarray_0vlt,temp2,'k','LineWidth',1.5)
plot(vltarray_0vlt,temp1+temp2,'b','LineWidth',1.5)

end


%title(['Lateral Deflection (mm), Front Pile'],'FontSize',11.0,'FontName','Cambria')
xlabel('Vlt (%)','FontSize',11.0,'FontName','Cambria')
ylabel('Load capacity (kN)','FontSize',11.0,'FontName','Cambria')
% set(gca,'YDir','reverse');
        xlim([0 5])
%%

out_Tshaft_0vlt=out_Tshaft_0vlt*10^3;
out_Qbase_0vlt=out_Qbase_0vlt*10^3;

flim2_vectortemp       =-1*out_Tshaft_0vlt;%Tensile         %maximum tensile force of plastic elements (<=0)
flim1_vectortemp       =+1*out_Tshaft_0vlt;%Compression Inf;%It needs to be set as a limit prssure rather than a force considering area corresponding to each node%2.5*10^5/(0.5*1.2);%       %maximum compressive force of plastic elements   (>=0) %     flim1       = c *Nc *sc +  0.5 gamma * B *Ng *sg;
flim1_vectortemp(end,:,:)=out_Qbase_0vlt+flim1_vectortemp(end,:,:);

flim2_vector_singlepile=flim2_vectortemp(:,1,1);
flim1_vector_singlepile=flim1_vectortemp(:,1,1);

flim2_vector_vlt=reshape(flim2_vectortemp,[],1,size(flim2_vectortemp,3));
flim2_vector_vlt=squeeze(flim2_vector_vlt);

flim1_vector_vlt=reshape(flim1_vectortemp,[],1,size(flim1_vectortemp,3));
flim1_vector_vlt=squeeze(flim1_vector_vlt);

% flim1_vector_vlt(1:nnodes:end,:)=flim1_vector_vlt(2:nnodes:end,:)/10;%to
% avoid having zero values / instead use coesion
% flim2_vector_vlt(1:nnodes:end,:)=flim2_vector_vlt(2:nnodes:end,:)/10;%to
% avoid having zero values / instead use coesion
%%


end



%%




% for ii=1:length(Vltp_vect)




for ii=1:length(Vltp_vect)

    if switch_soil==1
        [ Uffx(:,:,ii),Uffz(:,:,ii) ]=u_LON(z_global,x_global,ht,Rt,epsilon_vect(ii),0,0.4999999);

    elseif switch_soil==9000
    
    
    [ uxdeep2,uzdeep2 ]=u_verr_deep( z_global,x_global,ht,epsilon ,rho_deep,nis,switch_shape);
         Uffx(:,:,ii) =uxdeep2*(Vltp_vect(ii)/Vltp_vect(end));
         Uffz(:,:,ii) =uzdeep2*(Vltp_vect(ii)/Vltp_vect(end));

    elseif switch_soil==3
            Uffx=Sx;
            Uffz =Sz;
         Uffx(:,:,ii) =Sx*(Vltp_vect(ii)/Vltp_vect(end));
         Uffz(:,:,ii) =Sz*(Vltp_vect(ii)/Vltp_vect(end));
        

    elseif switch_soil==151
        
        
        [ Uffx(:,:,ii),Uffz(:,:,ii) ]=u_LON(z_global,x_global,ht,Rt,epsilon_vect(ii),0,0.4999999);
        
        tempu=Uffz(:,:,ii);
        tempu(tempu<0)=tempu(tempu<0)/5;
        Uffz(:,:,ii)=tempu;
        
    elseif switch_soil==220303
        load out_u_reg_coeff_cd20_id03c.mat
        load out_u_reg_coeffx_cd20_id03c.mat
        id=0.3;
        I_gcurve_vect = 0.5*ht*ones(size(Vltp_vect));
        [ Uffx(:,:,ii),Uffz(:,:,ii) ]=u_SE(z_global,x_global,ht,Rt,epsilon_vect(ii),id,reg_ca,reg_cb,reg_c1,reg_c2,reg_c3,reg_c4,reg_c5,reg_c6,reg_cax,reg_cbx,reg_c1x,reg_c2x,reg_c3x,reg_c4x,reg_c5x,reg_c6x)  ; 

    elseif switch_soil==22005
        load out_u_reg_coeff_cd20_id05.mat
        load out_u_reg_coeffx_cd20_id05.mat         
        
        id=0.9;
        I_gcurve_vect = 0.5*ht*ones(size(Vltp_vect));
        [ Uffx(:,:,ii),Uffz(:,:,ii) ]=u_SE(z_global,x_global,ht,Rt,epsilon_vect(ii),id,reg_ca,reg_cb,reg_c1,reg_c2,reg_c3,reg_c4,reg_c5,reg_c6,reg_cax,reg_cbx,reg_c1x,reg_c2x,reg_c3x,reg_c4x,reg_c5x,reg_c6x)  ; 
    elseif switch_soil==22009
        load out_u_reg_coeff_cd20_id09.mat
        load out_u_reg_coeffx_cd20_id09.mat         
        
        id=0.9;
        I_gcurve_vect = 0.5*ht*ones(size(Vltp_vect));
        [ Uffx(:,:,ii),Uffz(:,:,ii) ]=u_SE(z_global,x_global,ht,Rt,epsilon_vect(ii),id,reg_ca,reg_cb,reg_c1,reg_c2,reg_c3,reg_c4,reg_c5,reg_c6,reg_cax,reg_cbx,reg_c1x,reg_c2x,reg_c3x,reg_c4x,reg_c5x,reg_c6x)  ; 
        
    elseif switch_soil==24505
        load out_u_reg_coeff_cd45_id05.mat
        load out_u_reg_coeffx_cd45_id05.mat         
        
        id=0.9;
        I_gcurve_vect = 0.5*ht*ones(size(Vltp_vect));
        [ Uffx(:,:,ii),Uffz(:,:,ii) ]=u_SE(z_global,x_global,ht,Rt,epsilon_vect(ii),id,reg_ca,reg_cb,reg_c1,reg_c2,reg_c3,reg_c4,reg_c5,reg_c6,reg_cax,reg_cbx,reg_c1x,reg_c2x,reg_c3x,reg_c4x,reg_c5x,reg_c6x)  ; 
    
  
    
    
    
    elseif switch_soil==999

        Slp=0;%+.5;
        So=Soinput;%+.5;%0.1/8.74*20;

            
            
        deltaS=Slp-So;
         [ Uffx(:,:,ii),Uffz(:,:,ii) ]=u_KOR(z_global,x_global,So,deltaS,L,d,Vltp_vect(ii));
        
    end
        
    Uffx_vect(:,:,ii)=coeffux*Uffx(:,:,ii);
    Uffz_vect(:,:,ii)=coeffuz*Uffz(:,:,ii);
    Uffy_vect(:,:,ii)=zeros(nnodes,npali);

   
    Uffz(isnan(Uffz)) = 0 ;
	Uffx(isnan(Uffx)) = 0 ;
	Uffx_vect(isnan(Uffx_vect)) = 0 ;
	Uffz_vect(isnan(Uffz_vect)) = 0 ; 
    
end

%%
if exist('switch_negative_fric')==0
    switch_negative_fric=0;
end
if exist('settlement_history_surface')==0
    settlement_history_surface=0;
end


if switch_negative_fric==1


cont_temp=1;
settlement_history=zeros(nnodes,1);
settlement_history=settlement_history_surface*[1:-1/(nnodes-1):0]';  
for ii=1:(length(Vltp_vect)-1)/2
        cont_temp=cont_temp+1;
        Uffy_vect_temp(:,:,cont_temp)=zeros(nnodes,npali);
        Uffx_vect_temp(:,:,cont_temp)=zeros(nnodes,npali);
        Uffz_vect_temp(:,:,cont_temp)=repmat(settlement_history,1,npali)*(ii/((length(Vltp_vect)-1)/2));        
end


Uffy_vect_temp(:,:,(length(Vltp_vect)-1)/2+2:length(Vltp_vect))=Uffy_vect(:,:,3:2:end);
Uffx_vect_temp(:,:,(length(Vltp_vect)-1)/2+2:length(Vltp_vect))=Uffx_vect(:,:,3:2:end);
Uffz_vect_temp(:,:,(length(Vltp_vect)-1)/2+2:length(Vltp_vect))=Uffz_vect(:,:,3:2:end)+repmat(settlement_history,1,npali);

Uffy_vect=Uffy_vect_temp;
Uffz_vect=Uffz_vect_temp;
Uffx_vect=Uffx_vect_temp;

Sx=Uffx_vect(:,:,end);
Sy=Uffy_vect(:,:,end);
Sz=Uffz_vect(:,:,end);
end
%%

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



%% Pile stiffness degradation

if exist('switch_headdeg')==1

    betadeg=zeros(num_incr_TUNNEL,npali);
    
    degval=0.9;
%     
%     betadeg(1:num_incr_TUNNEL/2,1)=[0:degval/(num_incr_TUNNEL/2-1):degval];
%     betadeg(1:num_incr_TUNNEL/2,3)=[0:degval/(num_incr_TUNNEL/2-1):degval];
%     betadeg(1:num_incr_TUNNEL/2,4)=[0:degval/(num_incr_TUNNEL/2-1):degval];
% 
%     betadeg(num_incr_TUNNEL/2+1:end,1)=[degval];
%     betadeg(num_incr_TUNNEL/2+1:end,3)=[degval];
%     betadeg(num_incr_TUNNEL/2+1:end,4)=[degval];

    betadeg(1:num_incr_TUNNEL/2,1)=0;
    betadeg(1:num_incr_TUNNEL/2,3)=0;
    betadeg(1:num_incr_TUNNEL/2,4)=0;

    betadeg(num_incr_TUNNEL/2+1:end,1)=[0:degval/(num_incr_TUNNEL/2-1):degval];
    betadeg(num_incr_TUNNEL/2+1:end,3)=[0:degval/(num_incr_TUNNEL/2-1):degval];
    betadeg(num_incr_TUNNEL/2+1:end,4)=[0:degval/(num_incr_TUNNEL/2-1):degval];
    


    for ii=1:num_incr_TUNNEL


            betadeg_matrix=zeros(size(KKpalo));
            for jj=1:npali
            betadeg_matrix(:,1+6*nnodes*(jj-1):6*nnodes*jj)=betadeg(ii,jj);
            end
degscalar(1:50,1,1)=(1-betadeg(:,1));
% % %             SS_deg(:,:,ii)=+KKtors+KKpalo-betadeg_matrix.*KKpalo_negative;

    end
end





%% 
% =========================================================================
% ELASTO SOLUTION (KLAR ET AL 2005)
% =========================================================================


%Rename vectors
KKfoot=KKpalo;
nelementi_foot=nelementi;
nfoot=npali;
KKextra=KKtors;
KKframe_fixed_foot=KKframe_fixed;
KKframe_hinged_foot=KKframe_hinged;


% %GLOBAL STIFFNESS MATRIX    
%     KKpg=KKpalo+KKsoil+KKtors;  %KKsoil_base+
% %GLOBAL STIFFNESS MATRIX FRAME WITH FIXED PILE HEADS  
%     KKpg_fr_i=KKframe_fixed+KKpalo+KKsoil+KKtors;  %+KKsoil_base
% %GLOBAL STIFFNESS MATRIX FRAME WITH HINGED PILE HEADS  
%     KKpg_fr_h=KKframe_hinged+KKpalo+KKsoil+KKtors; %+KKsoil_base




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

lamdastd_hyp=lamdastd;  %initial value for the hyperbolic behaviour
Kst_hyp=Kst;            %initial value for the hyperbolic behaviour

% %%%%%%%%%%%%
% %debug
% % lamdasts=0*lamdasts;
% % % lamdastd(1:6:end,1:6:end)=lamdastd(1,1);
% % % lamdastd(2:6:end,2:6:end)=lamdastd(2,2);
% % % lamdastd(3:6:end,3:6:end)=lamdastd(3,3);
% % % lamdastd=diag(diag(lamdastd));
% %%%%%%%%%%%%5


%TOTAL STIFFNESS MATRIX
KK_el=KKfoot+Kst+KKextra+Kst*lamdasts*KKfoot;

KK_fr_i_el=(KKfoot+KKframe_fixed_foot)+Kst+KKextra+Kst*lamdasts*(KKfoot+KKframe_fixed_foot);
KK_fr_h_el=(KKfoot+KKframe_hinged_foot)+Kst+KKextra+Kst*lamdasts*(KKfoot+KKframe_hinged_foot);

% % if switch_cap==2
% % KK_fr_i_el=KK_fr_i_el+KKcap_fixed;
% % KK_fr_h_el=KK_fr_h_el+KKcap_fixed;
% % end


for jj=1:num_incr_TUNNEL
    
    %GF INPUT
    UCAT=Dfft_gruppo_vect(:,jj+1); %Tunnelling-induced greenfield movements
        
    %FORCE VECTOR
    %Tunnelling-induced forces
    F_el=Kst*UCAT;
    
    
    P_el=zeros(size(F_el));
    %External forces (concentrated at the centre of the footings)
%     P_el(1+(nelementi_foot/2*6):6*nnodes_foot:end,1)=Fx_foot_centr;
%     P_el(3+(nelementi_foot/2*6):6*nnodes_foot:end,1)=Fz_foot_centr;
%     P_el(5+(nelementi_foot/2*6):6*nnodes_foot:end,1)=My_foot_centr;

    Qten=sum(flim2_vector_singlepile);
    Qtot=sum(flim1_vector_singlepile);
    %External forces (pile head)
    
    
    
    P_el(1:6*nnodes_foot:end,1)=0;
    P_el(3:6*nnodes_foot:end,1)=Qtot/SF0;
    if exist('ecc_Ncap_x')
%     P_el(5,1)=npali*Qtot*ecc_Ncap_x;
%     P_el(4,1)=npali*Qtot*ecc_Ncap_y;
    P_el(5:6*nnodes_foot:end,1)=-Qtot/SF0*ecc_Ncap_x;
    P_el(4:6*nnodes_foot:end,1)=+Qtot/SF0*ecc_Ncap_y;
    else
    end
   
    
    if exist('alpha_Ncap_x')
%     P_el(5,1)=npali*Qtot*ecc_Ncap_x;
%     P_el(4,1)=npali*Qtot*ecc_Ncap_y;
    P_el(1:6*nnodes_foot:end,1)=Qtot/SF0*tan(alpha_Ncap_x/360*2*pi);
    P_el(2:6*nnodes_foot:end,1)=Qtot/SF0*tan(alpha_Ncap_y/360*2*pi);
    
    P_el(5:6*nnodes_foot:end,1)=-Qtot/SF0*ecc_Ncap_x-0.5*Qtot/SF0*tan(alpha_Ncap_x/360*2*pi);
    P_el(4:6*nnodes_foot:end,1)=+Qtot/SF0*ecc_Ncap_y+0.5*Qtot/SF0*tan(alpha_Ncap_y/360*2*pi);
    else
    end
   
    
    

    

    
    
    
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

    F_soil_fr_i_el_P_z(:,jj)=F_soil_fr_i_el_P(3:6:end,jj);
    F_soil_fr_i_el_P_x(:,jj)=F_soil_fr_i_el_P(1:6:end,jj);

    
    
    F_soil_fr_i_el_deltaT(:,jj)=F_soil_fr_i_el_T(:,jj)-F_soil_fr_i_el_P(:,jj);
	F_soil_fr_i_el_deltaT_z(:,jj)=F_soil_fr_i_el_deltaT(3:6:end,jj);
    F_soil_fr_i_el_deltaT_x(:,jj)=F_soil_fr_i_el_deltaT(1:6:end,jj);
    
    
end


vect_Utotf_fr_i_el_P=U_node_fr_i_el_P;
vect_Utotf_fr_i_el_T=U_node_fr_i_el_T;
vect_Utotp_fr_i_el_T=zeros(size(U_node_fr_i_el_P));
vect_Utotp_fr_i_el_P=zeros(size(U_node_fr_i_el_P));

vect_Utotf_free_el_P=U_node_free_el_P;
vect_Utotf_free_el_T=U_node_free_el_T;
vect_Utotp_free_el_T=zeros(size(U_node_free_el_P));
vect_Utotp_free_el_P=zeros(size(U_node_free_el_P));

for i=1:nnodes_foot;
    Uxpg_fr_i_el_P(i,1:1:nfoot)  =vect_Utotf_fr_i_el_P(1+(i-1)*6:nnodes_foot*6:end,end);
    Uypg_fr_i_el_P(i,1:1:nfoot)  =vect_Utotf_fr_i_el_P(2+(i-1)*6:nnodes_foot*6:end,end);
    Uzpg_fr_i_el_P(i,1:1:nfoot)  =vect_Utotf_fr_i_el_P(3+(i-1)*6:nnodes_foot*6:end,end);
    Ur2pg_fr_i_el_P(i,1:1:nfoot) =vect_Utotf_fr_i_el_P(5+(i-1)*6:nnodes_foot*6:end,end);
    
    Uxpg_fr_i_el_T(i,1:1:nfoot)  =vect_Utotf_fr_i_el_T(1+(i-1)*6:nnodes_foot*6:end,end);
    Uypg_fr_i_el_T(i,1:1:nfoot)  =vect_Utotf_fr_i_el_T(2+(i-1)*6:nnodes_foot*6:end,end);
    Uzpg_fr_i_el_T(i,1:1:nfoot)  =vect_Utotf_fr_i_el_T(3+(i-1)*6:nnodes_foot*6:end,end);
    Ur2pg_fr_i_el_T(i,1:1:nfoot) =vect_Utotf_fr_i_el_T(5+(i-1)*6:nnodes_foot*6:end,end);


    Uxpg_free_el_P(i,1:1:nfoot)  =vect_Utotf_free_el_P(1+(i-1)*6:nnodes_foot*6:end,end);
    Uypg_free_el_P(i,1:1:nfoot)  =vect_Utotf_free_el_P(2+(i-1)*6:nnodes_foot*6:end,end);
    Uzpg_free_el_P(i,1:1:nfoot)  =vect_Utotf_free_el_P(3+(i-1)*6:nnodes_foot*6:end,end);
    Ur2pg_free_el_P(i,1:1:nfoot) =vect_Utotf_free_el_P(5+(i-1)*6:nnodes_foot*6:end,end);
    
    Uxpg_free_el_T(i,1:1:nfoot)  =vect_Utotf_free_el_T(1+(i-1)*6:nnodes_foot*6:end,end);
    Uypg_free_el_T(i,1:1:nfoot)  =vect_Utotf_free_el_T(2+(i-1)*6:nnodes_foot*6:end,end);
    Uzpg_free_el_T(i,1:1:nfoot)  =vect_Utotf_free_el_T(3+(i-1)*6:nnodes_foot*6:end,end);
    Ur2pg_free_el_T(i,1:1:nfoot) =vect_Utotf_free_el_T(5+(i-1)*6:nnodes_foot*6:end,end);


end

Uxpg_fr_i_el_deltaT  =Uxpg_fr_i_el_T-Uxpg_fr_i_el_P;
Uypg_fr_i_el_deltaT  =Uypg_fr_i_el_T-Uypg_fr_i_el_P;
Uzpg_fr_i_el_deltaT  =Uzpg_fr_i_el_T-Uzpg_fr_i_el_P;
Ur2pg_fr_i_el_deltaT =Ur2pg_fr_i_el_T-Ur2pg_fr_i_el_P;

Uxpg_free_el_deltaT  =Uxpg_free_el_T-Uxpg_free_el_P;
Uypg_free_el_deltaT  =Uypg_free_el_T-Uypg_free_el_P;
Uzpg_free_el_deltaT  =Uzpg_free_el_T-Uzpg_free_el_P;
Ur2pg_free_el_deltaT =Ur2pg_free_el_T-Ur2pg_free_el_P;

vect_Utotp_fr_i_el_P_z=vect_Utotp_fr_i_el_P(3:6:end,:);
vect_Utotp_fr_i_el_P_x=vect_Utotp_fr_i_el_P(1:6:end,:);
vect_Utotf_fr_i_el_P_z=vect_Utotf_fr_i_el_P(3:6:end,:);
vect_Utotf_fr_i_el_P_x=vect_Utotf_fr_i_el_P(1:6:end,:);
vect_Utotp_fr_i_el_T_z=vect_Utotp_fr_i_el_T(3:6:end,:);
vect_Utotp_fr_i_el_T_x=vect_Utotp_fr_i_el_T(1:6:end,:);
vect_Utotf_fr_i_el_T_z=vect_Utotf_fr_i_el_T(3:6:end,:);
vect_Utotf_fr_i_el_T_x=vect_Utotf_fr_i_el_T(1:6:end,:);





vect_Utotp_fr_i_el_P_r2=vect_Utotp_fr_i_el_P(5:6:end,:);
vect_Utotf_fr_i_el_P_r2=vect_Utotf_fr_i_el_P(5:6:end,:);
vect_Utotp_fr_i_el_T_r2=vect_Utotp_fr_i_el_T(5:6:end,:);
vect_Utotf_fr_i_el_T_r2=vect_Utotf_fr_i_el_T(5:6:end,:);

vect_Utotf_fr_i_el_deltaT_z = vect_Utotf_fr_i_el_T_z - repmat(vect_Utotf_fr_i_el_P_z(:,end), [1, size(vect_Utotf_fr_i_el_T_z,2)]);
vect_Utotf_fr_i_el_deltaT_x = vect_Utotf_fr_i_el_T_x - repmat(vect_Utotf_fr_i_el_P_x(:,end), [1, size(vect_Utotf_fr_i_el_T_x,2)]);
vect_Utotf_fr_i_el_deltaT_r2 = vect_Utotf_fr_i_el_T_r2 - repmat(vect_Utotf_fr_i_el_P_r2(:,end), [1, size(vect_Utotf_fr_i_el_T_r2,2)]);

vect_Utotf_fr_i_el_deltaT = vect_Utotf_fr_i_el_T - repmat(vect_Utotf_fr_i_el_P(:,end), [1, size(vect_Utotf_fr_i_el_T,2)]);
vect_Utotf_free_el_deltaT = vect_Utotf_free_el_T - repmat(vect_Utotf_free_el_P(:,end), [1, size(vect_Utotf_free_el_T,2)]);


%%
% =========================================================================
% ELASTO-PLASTIC SOLUTION (KLAR ET AL 2007)
% =========================================================================

if switch_ep==0
vect_Utotf_fr_i_ep_P=vect_Utotf_fr_i_el_P;
vect_Utotf_fr_i_ep_T=vect_Utotf_fr_i_el_T;
vect_Utotf_free_ep_P=vect_Utotf_free_el_P;
vect_Utotf_free_ep_T=vect_Utotf_free_el_T;


vect_Utotf_fr_i_ep_deltaT_z=vect_Utotf_fr_i_el_deltaT_z;
vect_Utotf_fr_i_ep_deltaT_x=vect_Utotf_fr_i_el_deltaT_x;
vect_Utotf_fr_i_ep_deltaT_r2=vect_Utotf_fr_i_el_deltaT_r2;   

vect_Utotf_fr_i_ep_deltaT=vect_Utotf_fr_i_el_deltaT;
vect_Utotf_free_ep_deltaT=vect_Utotf_free_el_deltaT;

end



if switch_ep==1

% ------------------------------------------------------------------------- 
% FREE FOOTINGS
% -------------------------------------------------------------------------


SS=+KKextra+KKfoot;
SS_original=SS;
if exist('switch_headdeg')==1
    for ii=1:num_incr_TUNNEL
        SS_deg(:,:,ii)=++KKextra+KKfoot-betadeg_matrix.*KKpalo_negative;
    end
    SS=SS_deg(:,:,end);
end



% script_ep_solution
script_ep_solution_loadpath
% script_ep_solution_loadunloadP

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
% % if switch_cap==2
% % KKframe_fixed_foot=KKframe_fixed_foot+KKcap_fixed;
% % end
% SS=+KKextra+KKfoot+KKframe_fixed_foot;
SS=+KKextra+KKfoot+KKframe_fixed_foot;
SS_original=SS;
if exist('switch_headdeg')==1
    for ii=1:num_incr_TUNNEL
        SS_deg(:,:,ii)=+KKextra+KKfoot+KKframe_fixed_foot-betadeg_matrix.*KKpalo_negative;
    end
    SS=SS_deg(:,:,end);
end

% script_ep_solution
script_ep_solution_loadpath
% script_ep_solution_loadunloadP

vect_Utotf_fr_i_ep_P=vect_Utotf_P;
vect_Utotf_fr_i_ep_T=vect_Utotf_T;
vect_Utotp_fr_i_ep_T=cumsum(vect_Upinc_T,2);
vect_Utotp_fr_i_ep_P=cumsum(vect_Upinc_P,2);




for ii=1:num_incr_TUNNEL
    F_soil_fr_i_ep_T(:,ii)=SS*vect_Utotf_T(:,ii)-P_el;
    F_soil_fr_i_ep_P(:,ii)=SS*vect_Utotf_P(:,end)-P_el;
    F_soil_fr_i_ep_deltaT(:,ii)=F_soil_fr_i_ep_T(:,ii)-F_soil_fr_i_ep_P(:,ii);
end
F_soil_fr_i_ep_T_z=F_soil_fr_i_ep_T(3:6:end,:);
F_soil_fr_i_ep_T_x=F_soil_fr_i_ep_T(1:6:end,:);

F_soil_fr_i_ep_P_z=F_soil_fr_i_ep_P(3:6:end,:);
F_soil_fr_i_ep_P_x=F_soil_fr_i_ep_P(1:6:end,:);

F_soil_fr_i_ep_deltaT_z=F_soil_fr_i_ep_T_z-F_soil_fr_i_ep_P_z;
F_soil_fr_i_ep_deltaT_x=F_soil_fr_i_ep_T_x-F_soil_fr_i_ep_P_x;


% sigma_soil_fr_i_el_T_z=F_soil_fr_i_el_T_z./repmat(el_area,[1,num_incr_TUNNEL]);
% sigma_soil_fr_i_el_T_x=F_soil_fr_i_el_T_x./repmat(el_area,[1,num_incr_TUNNEL]);
% sigma_soil_fr_i_ep_T_z=F_soil_fr_i_ep_T_z./repmat(el_area,[1,num_incr_TUNNEL]);
% sigma_soil_fr_i_ep_T_x=F_soil_fr_i_ep_T_x./repmat(el_area,[1,num_incr_TUNNEL]);
% 
% sigma_soil_fr_i_el_deltaT_z=F_soil_fr_i_el_deltaT_z./repmat(el_area,[1,num_incr_TUNNEL]);
% sigma_soil_fr_i_el_deltaT_x=F_soil_fr_i_el_deltaT_x./repmat(el_area,[1,num_incr_TUNNEL]);




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

vect_Utotf_free_ep_deltaT   = vect_Utotf_free_ep_T   - repmat(vect_Utotf_free_ep_P  (:,end), [1, size(vect_Utotf_free_ep_T  ,2)]);
vect_Utotf_free_ep_deltaT_z = vect_Utotf_free_ep_T_z - repmat(vect_Utotf_free_ep_P_z(:,end), [1, size(vect_Utotf_free_ep_T_z,2)]);
vect_Utotf_free_ep_deltaT_x = vect_Utotf_free_ep_T_x - repmat(vect_Utotf_free_ep_P_x(:,end), [1, size(vect_Utotf_free_ep_T_x,2)]);
% vect_Utotf_free_ep_deltaT_r2 = vect_Utotf_free_ep_T_r2 - repmat(vect_Utotf_free_ep_P_r2(:,end), [1, size(vect_Utotf_free_ep_T_r2,2)]);

end



%%
% for i=1:nnodes;
% 	for ip=1:npali;
%   
%     Uxpg_fr_i_ep(i,ip)=vect_Utotf_fr_i_ep_T((ip-1)*(nnodes*6)+((i-1)*6)+1,end);
%     Uzpg_fr_i_ep(i,ip)=vect_Utotf_fr_i_ep_T((ip-1)*(nnodes*6)+((i-1)*6)+3,end);
% 
%     Uxpg_fr_i_ep_deltaT(i,ip)=vect_Utotf_fr_i_ep_deltaT((ip-1)*(nnodes*6)+((i-1)*6)+1,end);
%     Uzpg_fr_i_ep_deltaT(i,ip)=vect_Utotf_fr_i_ep_deltaT((ip-1)*(nnodes*6)+((i-1)*6)+3,end);
% 
%     Uxpg_free_ep(i,ip)=vect_Utotf_free_ep_T((ip-1)*(nnodes*6)+((i-1)*6)+1,end);
%     Uzpg_free_ep(i,ip)=vect_Utotf_free_ep_T((ip-1)*(nnodes*6)+((i-1)*6)+3,end);
% 
%     Uxpg_free_ep_deltaT(i,ip)=vect_Utotf_free_ep_deltaT((ip-1)*(nnodes*6)+((i-1)*6)+1,end);
%     Uzpg_free_ep_deltaT(i,ip)=vect_Utotf_free_ep_deltaT((ip-1)*(nnodes*6)+((i-1)*6)+3,end);
% 
%   
%     end
% end

for i=1:nnodes;
	for ip=1:npali;
  
    Uxpg_fr_i_ep_T(i,ip)=vect_Utotf_fr_i_ep_T((ip-1)*(nnodes*6)+((i-1)*6)+1,end);
    Uzpg_fr_i_ep_T(i,ip)=vect_Utotf_fr_i_ep_T((ip-1)*(nnodes*6)+((i-1)*6)+3,end);

    Uxpg_fr_i_ep_deltaT(i,ip)=vect_Utotf_fr_i_ep_deltaT((ip-1)*(nnodes*6)+((i-1)*6)+1,end);
    Uzpg_fr_i_ep_deltaT(i,ip)=vect_Utotf_fr_i_ep_deltaT((ip-1)*(nnodes*6)+((i-1)*6)+3,end);

    Uxpg_free_ep_T(i,ip)=vect_Utotf_free_ep_T((ip-1)*(nnodes*6)+((i-1)*6)+1,end);
    Uzpg_free_ep_T(i,ip)=vect_Utotf_free_ep_T((ip-1)*(nnodes*6)+((i-1)*6)+3,end);

    Uxpg_free_ep_deltaT(i,ip)=vect_Utotf_free_ep_deltaT((ip-1)*(nnodes*6)+((i-1)*6)+1,end);
    Uzpg_free_ep_deltaT(i,ip)=vect_Utotf_free_ep_deltaT((ip-1)*(nnodes*6)+((i-1)*6)+3,end);

  
    end
end

for i=1:nnodes;
	for ip=1:npali;
  
    Uxpg_fr_i_el_T(i,ip)=vect_Utotf_fr_i_el_T((ip-1)*(nnodes*6)+((i-1)*6)+1,end);
    Uzpg_fr_i_el_T(i,ip)=vect_Utotf_fr_i_el_T((ip-1)*(nnodes*6)+((i-1)*6)+3,end);

    Uxpg_fr_i_el_deltaT(i,ip)=vect_Utotf_fr_i_el_deltaT((ip-1)*(nnodes*6)+((i-1)*6)+1,end);
    Uzpg_fr_i_el_deltaT(i,ip)=vect_Utotf_fr_i_el_deltaT((ip-1)*(nnodes*6)+((i-1)*6)+3,end);

    Uxpg_free_el_T(i,ip)=vect_Utotf_free_el_T((ip-1)*(nnodes*6)+((i-1)*6)+1,end);
    Uzpg_free_el_T(i,ip)=vect_Utotf_free_el_T((ip-1)*(nnodes*6)+((i-1)*6)+3,end);

    Uxpg_free_el_deltaT(i,ip)=vect_Utotf_free_el_deltaT((ip-1)*(nnodes*6)+((i-1)*6)+1,end);
    Uzpg_free_el_deltaT(i,ip)=vect_Utotf_free_el_deltaT((ip-1)*(nnodes*6)+((i-1)*6)+3,end);

  
    end
end
% 
% for i=1:nnodes;
% 	for ip=1:npali;
%   
%     Uxpg_fr_i_el(i,ip)=vect_Utotf_fr_i_el_T((ip-1)*(nnodes*6)+((i-1)*6)+1,end);
%     Uzpg_fr_i_el(i,ip)=vect_Utotf_fr_i_el_T((ip-1)*(nnodes*6)+((i-1)*6)+3,end);
% 
%     Uxpg_fr_i_el_deltaT(i,ip)=vect_Utotf_fr_i_el_deltaT((ip-1)*(nnodes*6)+((i-1)*6)+1,end);
%     Uzpg_fr_i_el_deltaT(i,ip)=vect_Utotf_fr_i_el_deltaT((ip-1)*(nnodes*6)+((i-1)*6)+3,end);
% 
%     Uxpg_free_el(i,ip)=vect_Utotf_free_el_T((ip-1)*(nnodes*6)+((i-1)*6)+1,end);
%     Uzpg_free_el(i,ip)=vect_Utotf_free_el_T((ip-1)*(nnodes*6)+((i-1)*6)+3,end);
% 
%     Uxpg_free_el_deltaT(i,ip)=vect_Utotf_free_el_deltaT((ip-1)*(nnodes*6)+((i-1)*6)+1,end);
%     Uzpg_free_el_deltaT(i,ip)=vect_Utotf_free_el_deltaT((ip-1)*(nnodes*6)+((i-1)*6)+3,end);
% 
%   
%     end
% end


%%
%INNER FORCES
Ip=pi*(d/2)^4/4;

for i_p=1:npali
    for ii=1:nnodes
        if ii==1
            Mpg_fr_i_ep_T(ii,i_p)=KKframe_fixed_foot(5+6*(nnodes)*(i_p-1),:)*vect_Utotf_fr_i_ep_T(:,end);
            Mpg_free_ep_T(ii,i_p)=0;
            Mpg_fr_i_ep_deltaT(ii,i_p)=KKframe_fixed_foot(5+6*(nnodes)*(i_p-1),:)*vect_Utotf_fr_i_ep_deltaT(:,end);
            Mpg_free_ep_deltaT(ii,i_p)=0;
            
            
        elseif ii==nnodes

            Mpg_fr_i_ep_T(ii,i_p)=0;  
            Mpg_free_ep_T(ii,i_p)=0;
            Mpg_fr_i_ep_deltaT(ii,i_p)=0;  
            Mpg_free_ep_deltaT(ii,i_p)=0;
        else
            Mpg_fr_i_ep_T(ii,i_p)=       Ep*Ip*(Uxpg_fr_i_ep_T(ii-1,i_p)       -2*Uxpg_fr_i_ep_T(ii,i_p)       +Uxpg_fr_i_ep_T(ii+1,i_p))       /h_el^2;
            Mpg_free_ep_T(ii,i_p)=       Ep*Ip*(Uxpg_free_ep_T(ii-1,i_p)       -2*Uxpg_free_ep_T(ii,i_p)       +Uxpg_free_ep_T(ii+1,i_p))       /h_el^2;
            Mpg_fr_i_ep_deltaT(ii,i_p)=Ep*Ip*(Uxpg_fr_i_ep_deltaT(ii-1,i_p)-2*Uxpg_fr_i_ep_deltaT(ii,i_p)+Uxpg_fr_i_ep_deltaT(ii+1,i_p))/h_el^2;
            Mpg_free_ep_deltaT(ii,i_p)=Ep*Ip*(Uxpg_free_ep_deltaT(ii-1,i_p)-2*Uxpg_free_ep_deltaT(ii,i_p)+Uxpg_free_ep_deltaT(ii+1,i_p))/h_el^2;

            
 
        end
    end
end

for i_p=1:npali
    for ii=1:nnodes
        if ii==1
            Npg_fr_i_ep_T(ii,i_p)         =Ep*Ap*(Uzpg_fr_i_ep_T(ii+1,i_p)          -Uzpg_fr_i_ep_T(ii,i_p))       /h_el;
            Npg_free_ep_T(ii,i_p)         =Ep*Ap*(Uzpg_free_ep_T(ii+1,i_p)          -Uzpg_free_ep_T(ii,i_p))       /h_el;
            Npg_fr_i_ep_deltaT(ii,i_p)  =Ep*Ap*(Uzpg_fr_i_ep_deltaT(ii+1,i_p)   -Uzpg_fr_i_ep_deltaT(ii,i_p))/h_el;
            Npg_free_ep_deltaT(ii,i_p)  =Ep*Ap*(Uzpg_free_ep_deltaT(ii+1,i_p)   -Uzpg_free_ep_deltaT(ii,i_p))/h_el;

%             Npg_fr_i_ep_T(ii,i_p)=KKframe_fixed_foot(3+6*(nnodes)*(i_p-1),:)*vect_Utotf_fr_i_ep_T(:,end)-P_el(3:6*nnodes*npali:end);
%             Npg_free_ep_T(ii,i_p)=-P_el(3:6*nnodes*npali:end);
%             Npg_fr_i_ep_deltaT(ii,i_p)=KKframe_fixed_foot(3+6*(nnodes)*(i_p-1),:)*vect_Utotf_fr_i_ep_deltaT(:,end);
%             Npg_free_ep_deltaT(ii,i_p)=0;
%               
            
            
            if switch_ep==1
            Npg_fr_i_ep_P(ii,i_p)  =Ep*Ap*(Uzpg_fr_i_ep_P(ii+1,i_p)   -Uzpg_fr_i_ep_P(ii,i_p))/h_el;
            Npg_free_ep_P(ii,i_p)  =Ep*Ap*(Uzpg_free_ep_P(ii+1,i_p)   -Uzpg_free_ep_P(ii,i_p))/h_el;
%             
%             Npg_fr_i_ep_P(ii,i_p)=KKframe_fixed_foot(3+6*(nnodes)*(i_p-1),:)*vect_Utotf_fr_i_ep_P(:,end)-P_el(3:6*nnodes*npali:end);
%             Npg_free_ep_P(ii,i_p)=-P_el(3:6*nnodes*npali:end);
            
            
            end
            
        elseif ii==nnodes          
            Npg_fr_i_ep_T(ii,i_p)         =Ep*Ap*(Uzpg_fr_i_ep_T(ii,i_p)            -Uzpg_fr_i_ep_T(ii-1,i_p))       /h_el;    %note the Fztpg_spr_int
            Npg_free_ep_T(ii,i_p)         =Ep*Ap*(Uzpg_free_ep_T(ii,i_p)            -Uzpg_free_ep_T(ii-1,i_p))       /h_el;    %note the Fztpg_spr_int
            Npg_fr_i_ep_deltaT(ii,i_p)  =Ep*Ap*(Uzpg_fr_i_ep_deltaT(ii,i_p)     -Uzpg_fr_i_ep_deltaT(ii-1,i_p))/h_el;
            Npg_free_ep_deltaT(ii,i_p)  =Ep*Ap*(Uzpg_free_ep_deltaT(ii,i_p)     -Uzpg_free_ep_deltaT(ii-1,i_p))/h_el;
            if switch_ep==1
            Npg_fr_i_ep_P(ii,i_p)  =Ep*Ap*(Uzpg_fr_i_ep_P(ii,i_p)     -Uzpg_fr_i_ep_P(ii-1,i_p))/h_el;
            Npg_free_ep_P(ii,i_p)  =Ep*Ap*(Uzpg_free_ep_P(ii,i_p)     -Uzpg_free_ep_P(ii-1,i_p))/h_el;
            end
        else
            
            Npg_fr_i_ep_T(ii,i_p)         =Ep*Ap*(Uzpg_fr_i_ep_T(ii+1,i_p)          -Uzpg_fr_i_ep_T(ii-1,i_p))       /2/h_el;
            Npg_free_ep_T(ii,i_p)         =Ep*Ap*(Uzpg_free_ep_T(ii+1,i_p)          -Uzpg_free_ep_T(ii-1,i_p))       /2/h_el;
            Npg_fr_i_ep_deltaT(ii,i_p)  =Ep*Ap*(Uzpg_fr_i_ep_deltaT(ii+1,i_p)   -Uzpg_fr_i_ep_deltaT(ii-1,i_p))/2/h_el;
            Npg_free_ep_deltaT(ii,i_p)  =Ep*Ap*(Uzpg_free_ep_deltaT(ii+1,i_p)   -Uzpg_free_ep_deltaT(ii-1,i_p))/2/h_el;
            
             if switch_ep==1
            Npg_fr_i_ep_P(ii,i_p)  =     Ep*Ap*(Uzpg_fr_i_ep_P(ii+1,i_p)   -Uzpg_fr_i_ep_P(ii-1,i_p))/2/h_el;
            Npg_free_ep_P(ii,i_p)  =     Ep*Ap*(Uzpg_free_ep_P(ii+1,i_p)   -Uzpg_free_ep_P(ii-1,i_p))/2/h_el;
        
             end
        end
    end
end



for i_p=1:npali
    for ii=1:nnodes
        if ii==1
            Mpg_fr_i_el_T(ii,i_p)=KKframe_fixed_foot(5+6*(nnodes)*(i_p-1),:)*vect_Utotf_fr_i_el_T(:,end);
            Mpg_free_el_T(ii,i_p)=0;
            Mpg_fr_i_el_deltaT(ii,i_p)=KKframe_fixed_foot(5+6*(nnodes)*(i_p-1),:)*vect_Utotf_fr_i_el_deltaT(:,end);
            Mpg_free_el_deltaT(ii,i_p)=0;
            Mpg_free_el_deltaT_flex(ii,i_p)=0;
            
        elseif ii==nnodes

            Mpg_fr_i_el_T(ii,i_p)=0;  
            Mpg_free_el_T(ii,i_p)=0;
            Mpg_fr_i_el_deltaT(ii,i_p)=0;  
            Mpg_free_el_deltaT(ii,i_p)=0;
            Mpg_free_el_deltaT_flex(ii,i_p)=0;
        else
            Mpg_fr_i_el_T(ii,i_p)=       Ep*Ip*(Uxpg_fr_i_el_T(ii-1,i_p)       -2*Uxpg_fr_i_el_T(ii,i_p)       +Uxpg_fr_i_el_T(ii+1,i_p))       /h_el^2;
            Mpg_free_el_T(ii,i_p)=       Ep*Ip*(Uxpg_free_el_T(ii-1,i_p)       -2*Uxpg_free_el_T(ii,i_p)       +Uxpg_free_el_T(ii+1,i_p))       /h_el^2;
            Mpg_fr_i_el_deltaT(ii,i_p)=  Ep*Ip*(Uxpg_fr_i_el_deltaT(ii-1,i_p)-2*Uxpg_fr_i_el_deltaT(ii,i_p)+Uxpg_fr_i_el_deltaT(ii+1,i_p))/h_el^2;
            Mpg_free_el_deltaT(ii,i_p)=  Ep*Ip*(Uxpg_free_el_deltaT(ii-1,i_p)-2*Uxpg_free_el_deltaT(ii,i_p)+Uxpg_free_el_deltaT(ii+1,i_p))/h_el^2;
            Mpg_free_el_deltaT_flex(ii,i_p)=  Ep*Ip*(Sx(ii-1,i_p)-2*Sx(ii,i_p)+Sx(ii+1,i_p))/h_el^2;

 
        end
    end
end

for i_p=1:npali
    for ii=1:nnodes
        if ii==1
            Npg_fr_i_el_T(ii,i_p)         =Ep*Ap*(Uzpg_fr_i_el_T(ii+1,i_p)          -Uzpg_fr_i_el_T(ii,i_p))       /h_el;
            Npg_free_el_T(ii,i_p)         =Ep*Ap*(Uzpg_free_el_T(ii+1,i_p)          -Uzpg_free_el_T(ii,i_p))       /h_el;
% %             Npg_fr_i_el_T(ii,i_p)=KKframe_fixed_foot(3+6*(nnodes)*(i_p-1),:)*vect_Utotf_fr_i_el_T(:,end)-P_el(3:6*nnodes*npali:end);
% %             Npg_free_el_T(ii,i_p)=-P_el(3:6*nnodes*npali:end);
                     
            
            Npg_fr_i_el_deltaT(ii,i_p)  =Ep*Ap*(Uzpg_fr_i_el_deltaT(ii+1,i_p)   -Uzpg_fr_i_el_deltaT(ii,i_p))/h_el;
            Npg_free_el_deltaT(ii,i_p)  =Ep*Ap*(Uzpg_free_el_deltaT(ii+1,i_p)   -Uzpg_free_el_deltaT(ii,i_p))/h_el;
% %             Npg_fr_i_el_deltaT(ii,i_p)=KKframe_fixed_foot(3+6*(nnodes)*(i_p-1),:)*vect_Utotf_fr_i_el_deltaT(:,end);
% %             Npg_free_el_deltaT(ii,i_p)=0;
                        

            Npg_fr_i_el_P(ii,i_p)  =Ep*Ap*(Uzpg_fr_i_el_P(ii+1,i_p)   -Uzpg_fr_i_el_P(ii,i_p))/h_el;
            Npg_free_el_P(ii,i_p)  =Ep*Ap*(Uzpg_free_el_P(ii+1,i_p)   -Uzpg_free_el_P(ii,i_p))/h_el;
% %             Npg_fr_i_el_P(ii,i_p)=KKframe_fixed_foot(3+6*(nnodes)*(i_p-1),:)*vect_Utotf_fr_i_el_P(:,end)-P_el(3:6*nnodes*npali:end);
% %             Npg_free_el_P(ii,i_p)=-P_el(3:6*nnodes*npali:end);
                
            
        elseif ii==nnodes          
            Npg_fr_i_el_T(ii,i_p)         =Ep*Ap*(Uzpg_fr_i_el_T(ii,i_p)            -Uzpg_fr_i_el_T(ii-1,i_p))       /h_el;    %note the Fztpg_spr_int
            Npg_free_el_T(ii,i_p)         =Ep*Ap*(Uzpg_free_el_T(ii,i_p)            -Uzpg_free_el_T(ii-1,i_p))       /h_el;    %note the Fztpg_spr_int
            Npg_fr_i_el_deltaT(ii,i_p)  =Ep*Ap*(Uzpg_fr_i_el_deltaT(ii,i_p)     -Uzpg_fr_i_el_deltaT(ii-1,i_p))/h_el;
            Npg_free_el_deltaT(ii,i_p)  =Ep*Ap*(Uzpg_free_el_deltaT(ii,i_p)     -Uzpg_free_el_deltaT(ii-1,i_p))/h_el;

            Npg_fr_i_el_P(ii,i_p)  =Ep*Ap*(Uzpg_fr_i_el_P(ii,i_p)     -Uzpg_fr_i_el_P(ii-1,i_p))/h_el;
            Npg_free_el_P(ii,i_p)  =Ep*Ap*(Uzpg_free_el_P(ii,i_p)     -Uzpg_free_el_P(ii-1,i_p))/h_el;

        else
            
            Npg_fr_i_el_T(ii,i_p)         =Ep*Ap*(Uzpg_fr_i_el_T(ii+1,i_p)          -Uzpg_fr_i_el_T(ii-1,i_p))       /2/h_el;
            Npg_free_el_T(ii,i_p)         =Ep*Ap*(Uzpg_free_el_T(ii+1,i_p)          -Uzpg_free_el_T(ii-1,i_p))       /2/h_el;
            Npg_fr_i_el_deltaT(ii,i_p)  =Ep*Ap*(Uzpg_fr_i_el_deltaT(ii+1,i_p)   -Uzpg_fr_i_el_deltaT(ii-1,i_p))/2/h_el;
            Npg_free_el_deltaT(ii,i_p)  =Ep*Ap*(Uzpg_free_el_deltaT(ii+1,i_p)   -Uzpg_free_el_deltaT(ii-1,i_p))/2/h_el;
            

            Npg_fr_i_el_P(ii,i_p)  =     Ep*Ap*(Uzpg_fr_i_el_P(ii+1,i_p)   -Uzpg_fr_i_el_P(ii-1,i_p))/2/h_el;
            Npg_free_el_P(ii,i_p)  =     Ep*Ap*(Uzpg_free_el_P(ii+1,i_p)   -Uzpg_free_el_P(ii-1,i_p))/2/h_el;
        

        end
    end
end

%%

Area_shaft=[0.5*pi*d*h_el,pi*d*h_el*ones(1,nnodes-2),pi*d^2/4]';

Area_shaft_pg=repmat(Area_shaft, [1, npali]);


for i=1:nnodes
    for ip=1:npali
        index_direction=+3;     
        stresszpg_fr_i_el_P(i,ip)       =F_soil_fr_i_el_P((ip-1)*(nnodes*6)+((i-1)*6)+index_direction,end);
        stresszpg_fr_i_ep_P(i,ip)       =F_soil_fr_i_ep_P((ip-1)*(nnodes*6)+((i-1)*6)+index_direction,end);
       
        stresszpg_fr_i_el_T(i,ip)       =F_soil_fr_i_el_T((ip-1)*(nnodes*6)+((i-1)*6)+index_direction,end);
        stresszpg_fr_i_ep_T(i,ip)       =F_soil_fr_i_ep_T((ip-1)*(nnodes*6)+((i-1)*6)+index_direction,end);
        
        stresszpg_fr_i_el_deltaT(i,ip)  =F_soil_fr_i_el_deltaT((ip-1)*(nnodes*6)+((i-1)*6)+index_direction,end);
        stresszpg_fr_i_ep_deltaT(i,ip)  =F_soil_fr_i_ep_deltaT((ip-1)*(nnodes*6)+((i-1)*6)+index_direction,end);
        

        
    end
end

stresszpg_fr_i_el_P=-stresszpg_fr_i_el_P./Area_shaft_pg;
stresszpg_fr_i_ep_P=-stresszpg_fr_i_ep_P./Area_shaft_pg;
stresszpg_fr_i_el_T=-stresszpg_fr_i_el_T./Area_shaft_pg;
stresszpg_fr_i_ep_T=-stresszpg_fr_i_ep_T./Area_shaft_pg;
stresszpg_fr_i_el_deltaT=stresszpg_fr_i_el_deltaT./Area_shaft_pg;
stresszpg_fr_i_ep_deltaT=stresszpg_fr_i_ep_deltaT./Area_shaft_pg;


%%

Area_horizontal=[d*h_el,d*h_el*ones(1,nnodes-2),d*h_el]';

Area_horizontal_pg=repmat(Area_horizontal, [1, npali]);


for i=1:nnodes
    for ip=1:npali
        index_direction=+1;
        stressxpg_fr_i_el_P(i,ip)       =F_soil_fr_i_el_P((ip-1)*(nnodes*6)+((i-1)*6)+index_direction,end);
        stressxpg_fr_i_ep_P(i,ip)       =F_soil_fr_i_ep_P((ip-1)*(nnodes*6)+((i-1)*6)+index_direction,end);
       
        stressxpg_fr_i_el_T(i,ip)       =F_soil_fr_i_el_T((ip-1)*(nnodes*6)+((i-1)*6)+index_direction,end);
        stressxpg_fr_i_ep_T(i,ip)       =F_soil_fr_i_ep_T((ip-1)*(nnodes*6)+((i-1)*6)+index_direction,end);
        
        stressxpg_fr_i_el_deltaT(i,ip)  =F_soil_fr_i_el_deltaT((ip-1)*(nnodes*6)+((i-1)*6)+index_direction,end);
        stressxpg_fr_i_ep_deltaT(i,ip)  =F_soil_fr_i_ep_deltaT((ip-1)*(nnodes*6)+((i-1)*6)+index_direction,end);
        

        
    end
end

stressxpg_fr_i_el_P=stressxpg_fr_i_el_P./Area_horizontal_pg;
stressxpg_fr_i_ep_P=stressxpg_fr_i_ep_P./Area_horizontal_pg;
stressxpg_fr_i_el_T=stressxpg_fr_i_el_T./Area_horizontal_pg;
stressxpg_fr_i_ep_T=stressxpg_fr_i_ep_T./Area_horizontal_pg;
stressxpg_fr_i_el_deltaT=stressxpg_fr_i_el_deltaT./Area_horizontal_pg;
stressxpg_fr_i_ep_deltaT=stressxpg_fr_i_ep_deltaT./Area_horizontal_pg;
% 

%% Cap displacement

% Cap_result(iij,:)=[x_centre,x_centre/ht,U_cap_i_master,U_cap_i_master];
% U_fri_i_master_el=[mean(Uzpg_fr_i_el_deltaT(1,:)),mean(Uxpg_fr_i_el_deltaT(1,:)),mean(Ur2pg_fr_i_el_deltaT(1,:))];
% U_fri_i_master_ep=[mean(Uzpg_fr_i_ep_deltaT(1,:)),mean(Uxpg_fr_i_ep_deltaT(1,:)),mean(Ur2pg_fr_i_ep_deltaT(1,:)) ];
% Cap2_result(iij,:)=[x_centre,x_centre/ht,U_fri_i_master_el,U_fri_i_master_ep];

%% SIMPLIFIED METHODS FOR TUNNEL-PILE-STRUCTURE INTERACTION 
%

%Condensed stiffness matrix of structure
KK_fr_vspring=zeros(3*npali,3*npali);
if num_bay>0
    for i_p=1:npali
        KK_fr_vspring=RR_fixhead;
    end
end
%Condensed stiffness matrix of structure + vertical springs
%Interaction forces due to tunnelling
F_fnd_vspring=zeros(3*npali,1);
for i_p=1:npali
    KK_fr_vspring(3*(i_p-1)+2,3*(i_p-1)+2)=(nnodes-1)*h_el*kv+kb+KK_fr_vspring(3*(i_p-1)+2,3*(i_p-1)+2);

    for i=1:nnodes
        if i==nnodes
         F_fnd_vspring(3*(i_p-1)+2)=Sz(i,i_p)*(h_el/2*kv+kb)+F_fnd_vspring(3*(i_p-1)+2); 
        elseif i==1
         F_fnd_vspring(3*(i_p-1)+2)=Sz(i,i_p)*(h_el/2*kv+0)+F_fnd_vspring(3*(i_p-1)+2);
        else
         F_fnd_vspring(3*(i_p-1)+2)=Sz(i,i_p)*(h_el*kv+0)+F_fnd_vspring(3*(i_p-1)+2); 
        end
    end

end
%Displacement
%The first dof is constrainded (ux=0);
Upg_fr_vspring(1)=0;
Upg_fr_vspring(2:3*npali)=KK_fr_vspring(2:3*npali,2:3*npali)^(-1)*F_fnd_vspring(2:3*npali);
Uzpg_fr_vspring=Upg_fr_vspring(2:3:3*npali);
Uxpg_fr_vspring=Upg_fr_vspring(1:3:3*npali);
Urpg_fr_vspring=Upg_fr_vspring(3:3:3*npali);

%% CALCULATION OF THE REDUCTION FACTORS

%The following functions calculate relative deflection and horizontal
%strains of displacement profile at pile head

if exist('switch_modfc')==0
    switch_modfc=0;
end

if switch_modfc==0
else


% % %Relevant building dimensions
% % i_gf=Rt*(1.15*(ht/Rt)^0.9); %0.5*ht;
% % u_limit=2.5*i_gf;
% ind_limit=(abs(X))<u_limit;
ind_limit=(abs(0*X))==0;

% % %GF surface curve with 10cm space mesh
% % x_gfful=[min(X(ind_limit)):1:max(X(ind_limit))];
% % z_gfful=0*x_gfful;
% % [ Sxgfful,Szgfful ]=u_LON(z_gfful,x_gfful,ht,Rt,epsilon,0,0.49999);
% % Sxgfful=coeffux*Sxgfful;
% % Szgfful=coeffuz*Szgfful;

[ X_posinflpoint_gf,reldefl_mag_gf,reldefl_mag_sag_gf,reldefl_mag_hog_gf,l_sag_gf,l_hog_gf ] = f_reldeflection(0,X(ind_limit),Y(ind_limit),Sz(1,ind_limit));
% [ X_posinflpoint_gf,reldefl_mag_gf,reldefl_mag_sag_gf,reldefl_mag_hog_gf,l_sag_gf,l_hog_gf ] = f_reldeflection(1,x_gfful,z_gfful,Szgfful);
[ X_posinflpoint_free,reldefl_mag_free,reldefl_mag_sag_free,reldefl_mag_hog_free,l_sag_free,l_hog_free ] = f_reldeflection(0,X(ind_limit),Y(ind_limit),Uzpg_free_ep_deltaT(1,ind_limit));
[ X_posinflpoint_fr_i,reldefl_mag_fr_i,reldefl_mag_sag_fr_i,reldefl_mag_hog_fr_i,l_sag_fr_i,l_hog_fr_i ] = f_reldeflection(0,X(ind_limit),Y(ind_limit),Uzpg_fr_i_ep_deltaT(1,ind_limit));
% [ X_posinflpoint_fr_h,reldefl_mag_fr_h,reldefl_mag_sag_fr_h,reldefl_mag_hog_fr_h,l_sag_fr_h,l_hog_fr_h ] = f_reldeflection(1,X(ind_limit),Y(ind_limit),Uzpg_fr_i_ep_deltaT(1,ind_limit));
% % [ X_posinflpoint_fr_vspring,reldefl_mag_fr_vspring,reldefl_mag_sag_fr_vspring,reldefl_mag_hog_fr_vspring,l_sag_fr_vspring,l_hog_fr_vspring ] = f_reldeflection(1,X(ind_limit),Y(ind_limit),Uzpg_fr_vspring(1,ind_limit));

% X_posinflpoint=[ X_posinflpoint_gf; X_posinflpoint_free; X_posinflpoint_fr_i; X_posinflpoint_fr_h; X_posinflpoint_fr_vspring];

[ eps_ht_gf,eps_hc_gf,pos_ht_gf,pos_hc_gf] = f_hrzstr(0,X(ind_limit),Y(ind_limit),Sx(1,ind_limit)); %0-1-foundn
% [ eps_ht_gf,eps_hc_gf,pos_ht_gf,pos_hc_gf] = f_hrzstr(0,x_gfful,z_gfful,Sxgfful); %0-1-foundn
[ eps_ht_free,eps_hc_free,pos_ht_free,pos_hc_free] = f_hrzstr(0,X(ind_limit),Y(ind_limit),Uxpg_free_ep_deltaT(1,ind_limit));
[ eps_ht_fr_i,eps_hc_fr_i,pos_ht_fr_i,pos_hc_fr_i] = f_hrzstr(0,X(ind_limit),Y(ind_limit),Uxpg_fr_i_ep_deltaT(1,ind_limit));
% % [ eps_ht_fr_h,eps_hc_fr_h,pos_ht_fr_h,pos_hc_fr_h] = f_hrzstr(0,X(ind_limit),Y(ind_limit),Uxpg_fr_h(1,ind_limit));

%Without limit
% [ eps_ht_gf,eps_hc_gf,pos_ht_gf,pos_hc_gf] = f_hrzstr(0,X,Y,Sx(1,:)); %0-1-foundn
% [ eps_ht_free,eps_hc_free,pos_ht_free,pos_hc_free] = f_hrzstr(0,X,Y,Uxpg_free(1,:));
% [ eps_ht_fr_i,eps_hc_fr_i,pos_ht_fr_i,pos_hc_fr_i] = f_hrzstr(0,X,Y,Uxpg_fr_i(1,:));
% [ eps_ht_fr_h,eps_hc_fr_h,pos_ht_fr_h,pos_hc_fr_h] = f_hrzstr(0,X,Y,Uxpg_fr_h(1,:));


%Reduction factors for the relative deflection
red_factor_reldefl_hog_free=reldefl_mag_hog_free./reldefl_mag_hog_gf;
red_factor_reldefl_sag_free=reldefl_mag_sag_free./reldefl_mag_sag_gf;
red_factor_reldefl_hog_fr_i=reldefl_mag_hog_fr_i./reldefl_mag_hog_gf;
red_factor_reldefl_sag_fr_i=reldefl_mag_sag_fr_i./reldefl_mag_sag_gf;
% red_factor_reldefl_hog_fr_h=reldefl_mag_hog_fr_h./reldefl_mag_hog_gf;
% red_factor_reldefl_sag_fr_h=reldefl_mag_sag_fr_h./reldefl_mag_sag_gf;
% red_factor_reldefl_hog_fr_vspring=reldefl_mag_hog_fr_vspring./reldefl_mag_hog_gf;
% red_factor_reldefl_sag_fr_vspring=reldefl_mag_sag_fr_vspring./reldefl_mag_sag_gf;
red_factor_reldefl_hog_fr_h=red_factor_reldefl_hog_fr_i;
red_factor_reldefl_sag_fr_h=red_factor_reldefl_sag_fr_i;
red_factor_reldefl_hog_fr_vspring=red_factor_reldefl_hog_fr_i;
red_factor_reldefl_sag_fr_vspring=red_factor_reldefl_sag_fr_i;

%Reduction factors for the horizontal strains
red_factor_eht_free=eps_ht_free./eps_ht_gf;
red_factor_ehc_free=eps_hc_free./eps_hc_gf;
red_factor_eht_fr_i=eps_ht_fr_i./eps_ht_gf;
red_factor_ehc_fr_i=eps_hc_fr_i./eps_hc_gf;
% red_factor_eht_fr_h=eps_ht_fr_h./eps_ht_gf;
% red_factor_ehc_fr_h=eps_hc_fr_h./eps_hc_gf;
red_factor_ehc_fr_h=red_factor_ehc_fr_i;
red_factor_eht_fr_h=red_factor_eht_fr_i;

%Normalized and effective equivalent stiffness computed [Goh and Mair 
% (2014) equivalent coefficient is used to account for storey stiffness]
if foundn==1 && num_storey==0
    ro_sag= Eframe*(bf*df^3/12)/Es/l_sag_gf^3;   %it is slg*rosag
    ro_hog= Eframe*(bf*df^3/12)/Es/l_hog_gf^3;   %it is slg*rohog
    alpha_r =	Eframe*(bf*df)/Es/(max(X)-min(X));     
    
    EI_eq_fr_hog= Eframe*(bf*df^3/12);
    EI_eq_fr_sag= Eframe*(bf*df^3/12);
    EA_eq_fr=Eframe*(bf*df);
    
        
else
    kcfr=Eframe*(bc*dc^3/12)/spanz_frame;
    kbfr=Eframe*(bb*db^3/12)/spanx_frame;

    num_bay_hog=l_hog_gf/spanx_frame;
    num_bay_sag=l_sag_gf/spanx_frame;
    
    c_eq_fr_hog=1+num_bay_hog^2*(2*kcfr)/(2*kcfr+kbfr);
    c_eq_fr_sag=1+num_bay_sag^2*(2*kcfr)/(2*kcfr+kbfr);
    
    c_eq_fr_hog_topfloor=1+num_bay_hog^2*(kcfr)/(kcfr+kbfr);
    c_eq_fr_sag_topfloor=1+num_bay_sag^2*(kcfr)/(kcfr+kbfr);
    
% 	c_eq_fr_hog=1;
%     c_eq_fr_sag=1;
%     
%     c_eq_fr_hog_topfloor=1;
%     c_eq_fr_sag_topfloor=1;

    % %Equivalent square beam at foundation level modelling frame
    if foundn==1
        EI_eq_fr_hog= Eframe*(bf*df^3/12+(num_storey-1)*c_eq_fr_hog*bb*db^3/12+c_eq_fr_hog_topfloor*bb*db^3/12);
        EI_eq_fr_sag= Eframe*(bf*df^3/12+(num_storey-1)*c_eq_fr_sag*bb*db^3/12+c_eq_fr_sag_topfloor*bb*db^3/12);
        EA_eq_fr=Eframe*(bf*df+num_storey*bb*db);
        alpha_r =	Eframe*(bf*df)/Es/spanx_frame;   
        alpha_r = alpha_r + 1/Es*(3*kbfr*kcfr/spanz_frame^2/(2*kbfr+3*kcfr)); 
    else
        EI_eq_fr_hog= Eframe*((num_storey-1)*c_eq_fr_hog*bb*db^3/12+c_eq_fr_hog_topfloor*bb*db^3/12);
        EI_eq_fr_sag= Eframe*((num_storey-1)*c_eq_fr_sag*bb*db^3/12+c_eq_fr_sag_topfloor*bb*db^3/12);
        EA_eq_fr=Eframe*(num_storey*bb*db);
        alpha_r =	1/Es*(3*kbfr*kcfr/spanz_frame^2/(2*kbfr+3*kcfr));    
    end
    ro_sag=EI_eq_fr_sag/Es/l_sag_gf^3;
    ro_hog=EI_eq_fr_hog/Es/l_hog_gf^3;
    
end

result_red_factor=[...
	Es,...
    Eb,...
    nis,...
    nib,...
    Ep,...
    Eframe,...
    num_storey,...%7
    npali,...
    spanx_frame,...
    x_centre,...
    d,...
    L, ...
    ht,...
    Rt,...%14
	EI_eq_fr_sag,...
    EI_eq_fr_hog,...
    EA_eq_fr,...
    ro_sag,... %18
    ro_hog, ...
    alpha_r,...
    reldefl_mag_sag_gf,...
    reldefl_mag_hog_gf,...
    red_factor_reldefl_sag_free,...%23
    red_factor_reldefl_hog_free,...
    red_factor_reldefl_sag_fr_i,...%25
    red_factor_reldefl_hog_fr_i,...
    red_factor_reldefl_sag_fr_h,...
    red_factor_reldefl_hog_fr_h,...
    red_factor_reldefl_sag_fr_vspring,...
    red_factor_reldefl_hog_fr_vspring,...
    red_factor_eht_free,...
    red_factor_ehc_free,...
    red_factor_eht_fr_h,...
    red_factor_ehc_fr_h,...
    red_factor_eht_fr_i,...%35
    red_factor_ehc_fr_i,...
    l_sag_gf,...
    l_hog_gf,...
    eps_ht_gf,...%39
    eps_hc_gf];


%% SAVE
if exist('Data_red_factor.mat')==0 
	data_red_factor(1,:)=result_red_factor;
	save('Data_red_factor.mat','data_red_factor')  
else
    load Data_red_factor.mat
	i=1+size(data_red_factor,1);
	data_red_factor(i,:)=result_red_factor;
    save('Data_red_factor.mat','data_red_factor')  
end
end