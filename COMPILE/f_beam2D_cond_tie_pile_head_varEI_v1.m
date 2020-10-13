function [ RR ] = f_beam2D_cond_tie_pile_head_varEI_v1( xx,yy,EI_mat,EA_mat)
% This code analyse 2D frame. All units are in N,m
% counter-clockwise FEM is +m and rot
%right and top is +u anf f


%% Inputs

GDof=3*length(xx);
Elements_number=length(xx)-1;
Bays_Number=length(xx)-1;
%% 

%Global stiffnes matrix
KG=zeros(GDof);
for i=1:Elements_number
    EA=EA_mat(i);
    EI=EI_mat(i);
    
    %Calculation of element length & angle of inclination
    x1=xx(i);
    y1=yy(i);
    x2=xx(i+1);
    y2=yy(i+1);
    
    le=sqrt((x2-x1)^2); %Length of the memeber
    Le(i)=le;

    
    %Element stiffness
    Ke=[EA/le  ,0          ,0          ,-EA/le ,0          ,0;
        0      ,12*EI/le^3 ,6*EI/le^2  ,0      ,-12*EI/le^3 ,6*EI/le^2;
        0      ,6*EI/le^2  ,4*EI/le    ,0      ,-6*EI/le^2 ,2*EI/le;
        -EA/le ,0          ,0          ,EA/le  ,0          ,0;
        0      ,-12*EI/le^3 ,-6*EI/le^2 ,0      ,12*EI/le^3 ,-6*EI/le^2;
        0      ,6*EI/le^2  ,2*EI/le    ,0      ,-6*EI/le^2 ,4*EI/le];
    


    %Assembly of KG
    KG(1+3*(i-1):6+3*(i-1),1+3*(i-1):6+3*(i-1))= KG(1+3*(i-1):6+3*(i-1),1+3*(i-1):6+3*(i-1))+Ke;
end
        
        
%% 


% Number of piles/supports
num_supp=Bays_Number+1;


   
for iip=1:num_supp
for iidof=1:3

    [  Reactions ] =KG(:,3*(iip-1)+iidof);


    %RReactions condensed stiffness matrix written according to 
    %- displacement constraint convention of the tunnel-pile interaction method (pos if
    % right down counter)
    %- structural forces convention of structural analyses (pos right top counter)

    if iidof==1
        RReactions(:,num_supp*(2-1)+iip)=Reactions; %Matrix reaction - Stiffness matrix
    end
    if iidof==2
        RReactions(:,num_supp*(1-1)+iip)=-1*Reactions; %Matrix reaction - Stiffness matrix (different sign convention)
    end
    if iidof==3
        RReactions(:,num_supp*(3-1)+iip)=Reactions; %Matrix reaction - Stiffness matrix
    end



end
end

    %RR condensed stiffness matrix written according to the convention of the
    %tunnel-pile interaction method
    %- displacement constraint convention of the tunnel-pile interaction method (pos if
    % right down counter)
    %- structural forces convention of structural analyses (pos left down counter)
    for iip=1:num_supp
    for iidof=1:3
        if iidof==1
        RR(3*(iip-1)+2,:)=-1*RReactions(3*(iip-1)+iidof,:); % Rearrenged degree of freedom
        end

        if iidof==2
        RR(3*(iip-1)+1,:)=-1*RReactions(3*(iip-1)+iidof,:); % Rearrenged degree of freedom
        end

        if iidof==3
        RR(3*(iip-1)+3,:)=RReactions(3*(iip-1)+iidof,:); % Rearrenged degree of freedom
        end

    end
    end

end

