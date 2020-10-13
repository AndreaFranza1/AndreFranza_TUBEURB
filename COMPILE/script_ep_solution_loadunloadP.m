UINC_iprevious=zeros(6*nfoot*(nelementi_foot+1),1);
UINCIP_iprevious=zeros(6*nfoot*(nelementi_foot+1),1);
% % dUCAT2=Dfft_gruppo/num_incr_TUNNEL;
dUCAT_vect=diff(Dfft_gruppo_vect,1,2);
dUCAT_vect_memory=dUCAT_vect;

%% Solution for External load and selfweight

if sum(P_el~=zeros(size(P_el)))~=0
    i_inc=0;
    num_incr=num_incr_P;
    
    UINC_iprevious=zeros(6*nfoot*(nelementi_foot+1),1);
    P_elINC_iprevious=zeros(6*nfoot*(nelementi_foot+1),1);
    
    dUCAT_vect=0.*dUCAT_vect;
    dP_el=P_el/num_incr_P;

%     script_ep_iteration2_2_mod_loadunloadP %this code improved by considering nedative df
%     

    if switch_nep==1
        script_ep_iteration2_2_mod_vect_hyp_unloadingP %this consider hyperbolic E degradation in z direction
    else
        script_ep_iteration2_2_mod_vect_unloadingP
    end 
    
    vect_Utotf_P=vect_UINC;
    vect_Upinc_P=vect_dUIP;
else
   	vect_Utotf_P=zeros(6*nfoot*(nelementi_foot+1),1);
    vect_Upinc_P=zeros(6*nfoot*(nelementi_foot+1),1);

end

%% Solution for Tunnelling-induced effects
if Vltp~=0
    i_inc=0;
    num_incr=num_incr_TUNNEL;
    
	UINC_iprevious=vect_Utotf_P(:,end);
	P_elINC_iprevious=0*P_el(:,end);
% % % %     P_elINC_iprevious=P_el(:,end);
    
    dUCAT_vect=dUCAT_vect_memory;
	dP_el=zeros(6*nfoot*(nelementi_foot+1),1);

%     script_ep_iteration2_2
    %     script_ep_iteration2_2_mod %this code improved by considering nedative df and unloading path

    if switch_nep==1
        script_ep_iteration2_2_mod_vect_hyp %this consider hyperbolic E degradation in z direction
    else
        script_ep_iteration2_2_mod_vect %this consider base and shaft capacity= f(z)
    end    

    vect_Utotf_T=vect_UINC;
    vect_Upinc_T=vect_dUIP;
end
