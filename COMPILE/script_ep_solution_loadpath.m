if exist('switch_sliderx')==0
switch_sliderx=999;
end



UINC_iprevious=zeros(6*nfoot*(nelementi_foot+1),1);
UINCIP_iprevious=zeros(6*nfoot*(nelementi_foot+1),1);
% % dUCAT2=Dfft_gruppo/num_incr_TUNNEL;
dUCAT_vect=diff(Dfft_gruppo_vect,1,2);
dUCAT_vect_memory=dUCAT_vect;

flim1_maxhistory=zeros(1*nfoot*(nelementi_foot+1),1);
flim2_maxhistory=zeros(1*nfoot*(nelementi_foot+1),1);
%% Solution for External load and selfweight

if sum(P_el~=zeros(size(P_el)))~=0
    i_inc=0;
    num_incr=num_incr_P;
    
    UINC_iprevious=zeros(6*nfoot*(nelementi_foot+1),1);
    P_elINC_iprevious=zeros(6*nfoot*(nelementi_foot+1),1);
    
    dUCAT_vect=0.*dUCAT_vect;
    

    %%%%%%%%%%%
    dP_el=P_el/num_incr_P;
    if mod(num_incr,4)==1   % odd number
        stop
    end
    if switch_displfound==1  %displacement pile
        Percentage_preloading=0.95;%of the total load Qtot
        dP_el_vect=[(P_el*SF0)*Percentage_preloading/num_incr_P*2*ones(size(P_el,2),num_incr/2),....
                   -(P_el*SF0)*Percentage_preloading/num_incr_P*4*ones(size(P_el,2),num_incr/4),...
                                             +(P_el)/num_incr_P*4*ones(size(P_el,2),num_incr/4)];
        P_el_vect=cumsum(dP_el_vect,2);
    elseif switch_displfound==2  %non displacement pile
        dP_el_vect=P_el./num_incr_P*ones(size(P_el,2),num_incr);
        P_el_vect=cumsum(dP_el_vect,2);

    end
    %%%%%%%%%%%%%%

    
%     script_ep_iteration2_2_mod %this code improved by considering nedative df and unloading path

    if switch_nep==1
%         script_ep_iteration2_2_mod_vect_hypv2 %this consider hyperbolic E degradation in z direction
        if switch_sliderx==1
        script_ep_iteration2_2_mod_vect_hypv2_limxz %this consider base and shaft capacity= f(z), also there is horizontal slider
        else
        script_ep_iteration2_2_mod_vect_hypv2%this consider base and shaft capacity= f(z)
        end  
        
    else
        if switch_sliderx==1
        script_ep_iteration2_2_mod_vect_limzx%this consider base and shaft capacity= f(z), also there is horizontal slider
        else
        script_ep_iteration2_2_mod_vect%this consider base and shaft capacity= f(z)
        end
%         script_ep_iteration2_2_mod_vect%this consider base and shaft capacity= f(z)
    end

    vect_Utotf_P=vect_UINC;
    vect_Upinc_P=vect_dUIP;
else
   	vect_Utotf_P=zeros(6*nfoot*(nelementi_foot+1),1);
    vect_Upinc_P=zeros(6*nfoot*(nelementi_foot+1),1);

end
if switch_nep==1        
red_fac_P=red_fac;
end
%% Activating the vertical spring

if switch_headspring==1
    
    KKheadspring=zeros(size(KKfoot));

    for ij=1:npali
    KKheadspring(3+(ij-1)*(6*nnodes),3+(ij-1)*(6*nnodes))=Kpile_el;
    end

    SS=SS+KKheadspring;    
end

%% Solution for Tunnelling-induced effects
if Vltp~=0
    i_inc=0;
    num_incr=num_incr_TUNNEL;
    
	UINC_iprevious=vect_Utotf_P(:,end);
	P_elINC_iprevious=P_el(:,end);
    
    dUCAT_vect=dUCAT_vect_memory;
	

    %%%%%%%%%%%
    dP_el=zeros(6*nfoot*(nelementi_foot+1),1);
    if switch_displfound==1  %displacement pile
        dP_el_vect=zeros(6*nfoot*(nelementi_foot+1),num_incr);
    elseif switch_displfound==2  %non displacement pile
        dP_el_vect=zeros(6*nfoot*(nelementi_foot+1),num_incr);
            elseif switch_displfound==3  %non displacement pile
        dP_el_vect=zeros(6*nfoot*(nelementi_foot+1),num_incr);
    end
    %%%%%%%%%%%%%%
    
    
    
    
%     script_ep_iteration2_2_mod %this code improved by considering nedative df and unloading path
    if switch_nep==1
%         script_ep_iteration2_2_mod_vect_hypv2 %this consider hyperbolic E degradation in z direction 
        if switch_sliderx==1
        script_ep_iteration2_2_mod_vect_hypv2_limxz%this consider base and shaft capacity= f(z), also there is horizontal slider
        else
        script_ep_iteration2_2_mod_vect_hypv2%this consider base and shaft capacity= f(z)
        end          
    else
        
        if switch_sliderx==1
        script_ep_iteration2_2_mod_vect_limzx%this consider base and shaft capacity= f(z), also there is horizontal slider
        else
        script_ep_iteration2_2_mod_vect%this consider base and shaft capacity= f(z)
        end
%         script_ep_iteration2_2_mod_vect%this consider base and shaft capacity= f(z)     
    end
    
    
    vect_Utotf_T=vect_UINC;
    vect_Upinc_T=vect_dUIP;
end
