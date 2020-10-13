function [ RReactions ] = f_frame2D_cond_tie_pile_head_v1( Bays_Number,Storeys_Number,lx,ly,E,bc,dc,bb,db,bf,df,foundn )
% This code analyse 2D frame. All units are in N,m
% counter-clockwise FEM is +m and rot
%right and top is +u anf f


%% Inputs

% Uniform distributed load
wb=0; %UDL on beams N/m
% Number of piles/supports
num_supp=Bays_Number+1;




%support=[support_Node No. ux uy uz; ....] & (0: restrained; 1: unrestrained)
% supports=[1 0 0 0;
%           2 0 0 0;
%           3 0 0 0]; %write the nodes number @ supports
for ip=1:num_supp
	supports(ip,:)=[ip, zeros(1,3)];
end

      

   
for iip=1:num_supp
for iidof=1:3

    % Settlment =[support_Node No. ux uy uz; ....]
    for ip=1:num_supp
        support_settlement(ip,:)=[ip, zeros(1,3)];
        % support_settlement=[1 0 1 0;
        %                     2 0 0 0;
        %                     3 0 0 0]; %write the nodes number @ supports settlement
    end
    support_settlement(iip,iidof+1)=1;

    [  Reactions ] = f_frame_script( iip,iidof,Bays_Number,Storeys_Number,lx,ly,E,bc,dc,bb,db,bf,df,support_settlement,supports,foundn,wb );








    %RReactions condensed stiffness matrix written according to 
    %- displacement constraint convention of the tunnel-pile interaction method (pos if
    % right down counter)
    %- structural forces convention of structural analyses (pos right top counter)
    if iidof==1
            sign_conve(1:3:num_supp*3,1)=+1;
            sign_conve(2:3:num_supp*3,1)=-1; 
            sign_conve(3:3:num_supp*3,1)=+1; 
        RReactions(:,(iip-1)*3+iidof)=Reactions.*sign_conve; %Matrix reaction - Stiffness matrix
    end
    if iidof==2
            sign_conve(1:3:num_supp*3,1)=-1;
            sign_conve(2:3:num_supp*3,1)=+1; 
            sign_conve(3:3:num_supp*3,1)=-1; 
        RReactions(:,(iip-1)*3+iidof)=Reactions.*sign_conve; %Matrix reaction - Stiffness matrix
    end
    if iidof==3
            sign_conve(1:3:num_supp*3,1)=+1;
            sign_conve(2:3:num_supp*3,1)=-1; 
            sign_conve(3:3:num_supp*3,1)=+1; 
        RReactions(:,(iip-1)*3+iidof)=Reactions.*sign_conve; %Matrix reaction - Stiffness matrix
    end

     




    end
end

    


%     %RR condensed stiffness matrix written according to the convention of the
%     %tunnel-pile interaction method
%     %- displacement constraint convention of the tunnel-pile interaction method (pos if
%     % right down counter)
%     %- structural forces convention of structural analyses (pos left down counter)
%     for i=1:num_supp
%      for j=1:num_supp   
%     KKframe(:,iip)=RReactions(:,iip)
%     KKframe=RReactions(:,iip+3*num_supp)
%     KKframe=RReactions(:,iip+6*num_supp)
%     
%      end
%     end

end

