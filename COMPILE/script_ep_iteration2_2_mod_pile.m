%This script considers both the plasticity of springs in direction 1 and 3
%(horizontal and vertical directions for the footings)

clear vect_UINC vect_dUIP

for ii=1:num_incr
    %Incremental step
    dUIP=zeros(6*nfoot*(nelementi_foot+1),1);
    dUCAP=zeros(6*nfoot*(nelementi_foot+1),1);
    dUCAT=dUCAT_vect(:,ii);
    while 1 %This while loop solve dUIP for a given dUCAP and dU
        
        dU=(SS+Kst)\(dP_el+Kst*dUCAP+Kst*dUCAT+Kst*dUIP);
        
        UINC=UINC_iprevious+dU;
        P_elINC=P_elINC_iprevious+dP_el;
        
        dF=-SS*dU+dP_el;
            %Plasticity conditions
              
% %             react_inc=-SS*UINC+P_elINC;
% %             ind_ip_3_up=react_inc(3:6:end)<flim2;%if f=-SS*UINC+P_elINC> 0 spring is compressed
% %             ind_ip_3_dw=react_inc(3:6:end)>flim1;
% %             ind_ip_3= ind_ip_3_up|ind_ip_3_dw ;
% %             dF(3:6:end)=dF(3:6:end).*(1-ind_ip_3);
% %             ind_ip_1=abs(react_inc(1:6:end))>mu*react_inc(3:6:end);
% %             dF(1:6:end)=dF(1:6:end).*(1-ind_ip_1);


            react_iprevious=    -SS*UINC_iprevious  +P_elINC_iprevious;
            react_inc=          -SS*UINC            +P_elINC;
            
            react_iprevious_3dof=react_iprevious(3:6:end);
            react_inc_3dof=react_inc(3:6:end);
            react_iprevious_1dof=react_iprevious(1:6:end);
            react_inc_1dof=react_inc(1:6:end);
            
            ind_ip_3_up=react_inc_3dof<flim2;%if f=-SS*UINC+P_elINC> 0 spring is compressed
            ind_ip_3_dw=react_inc_3dof>flim1;
            
            dF_3dof=dF(3:6:end);
            dF_1dof=dF(1:6:end);
            
            dF_3dof(ind_ip_3_up)=flim2-react_iprevious_3dof(ind_ip_3_up);
            dF_3dof(ind_ip_3_dw)=flim1-react_iprevious_3dof(ind_ip_3_dw); 
            dF(3:6:end)=dF_3dof;

%                    if max(ind_ip_3_up)>0
%                 pause
%             end     
            
            
            if ii==0
            fmax_1dof=mu*max(react_inc_3dof,0);
            ind_ip_1_up=react_inc_1dof> fmax_1dof ;
            ind_ip_1_dw=react_inc_1dof<-fmax_1dof ;
                        
            dF_1dof(ind_ip_1_up)= 0;
            dF_1dof(ind_ip_1_dw)= 0;         
            dF(1:6:end)=dF_1dof;

            else
            fmax_1dof=mu*max(react_inc_3dof,0);
            ind_ip_1_up=react_inc_1dof> fmax_1dof ;
            ind_ip_1_dw=react_inc_1dof<-fmax_1dof ;
                        
            dF_1dof(ind_ip_1_up)= fmax_1dof(ind_ip_1_up)-react_iprevious_1dof(ind_ip_1_up);
            dF_1dof(ind_ip_1_dw)=-fmax_1dof(ind_ip_1_dw)-react_iprevious_1dof(ind_ip_1_dw);         
            dF(1:6:end)=dF_1dof;
            end
            
        CI=((lamdasts+lamdastd)*dF+dUCAT+dUIP)./dU;
        CI(isnan(CI(:,1)),:)=[1];CI(2:6:end,:)=[1];CI(4:6:end,:)=[1];CI(5:6:end,:)=[1];CI(6:6:end,:)=[1];
        
        dUIP_old=dUIP;
        
        dUd=(lamdastd)*dF;
        dUIP(3:6:end)=dU(3:6:end)-dUd(3:6:end)-dUCAP(3:6:end)-dUCAT(3:6:end);
        dUIP(1:6:end)=dU(1:6:end)-dUd(1:6:end)-dUCAP(1:6:end)-dUCAT(1:6:end);
        
        dUCAP=((lamdasts)*dF+beta*dUCAP)/(1+beta);
               
        perr=(dUIP-dUIP_old)./(dUIP);
        perr(isnan(perr(:,1)),:)=[0];
        
        if max(abs(CI-1))<delta_err && max(abs(perr))<eps_err
            i_inc=i_inc+1;
            vect_UINC(:,i_inc)=UINC;
            vect_dUIP(:,i_inc)=dUIP;
            UINC_iprevious=UINC;
            P_elINC_iprevious=P_elINC;
            break
        end
        
    end
end