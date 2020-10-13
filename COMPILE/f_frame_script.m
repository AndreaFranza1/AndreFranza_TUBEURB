function [ Reactions ] = f_frame_script( iip,iidof,Bays_Number,Storeys_Number,lx,ly,E,bc,dc,bb,db,bf,df,support_settlement,supports,foundn,wb )
                
%% =======================================================================
% Geomety formation
if foundn==1;
    foundation_member=length(supports(:,1))-1;
else
    foundation_member=0;
end
Nodes_Number=(Bays_Number+1)*(Storeys_Number+1);
columns_Number=(Bays_Number+1)*Storeys_Number;
Beams_Number=Bays_Number*Storeys_Number;
Elements_Number=columns_Number+Beams_Number+foundation_member;
EAc=E*(bc*dc); EIc=E*(bc*dc^3/12);      % Columns 
EAb=E*(bb*db); EIb=E*(bb*db^3/12);      % Beams 
EAf=E*(bf*df); EIf=E*(bf*df^3/12);      % foundation beams

%Nodes_cordinates
counterN=0;
for ss=0:Storeys_Number
    for bb=0:Bays_Number
        counterN=counterN+1;
        Nodes_cordinates1(counterN,1)=counterN;
        Nodes_cordinates1(counterN,2)=lx*bb;
        Nodes_cordinates1(counterN,3)=ly*ss;
    end
end

Nodes_cordinates=Nodes_cordinates1(:,2:3);
%Connectivity matrix: columns are first and beams follows
Elements(:,1)=1:Elements_Number;
%colums
Elements(1:columns_Number,2)=Nodes_cordinates1(1:Storeys_Number*(Bays_Number+1),1);
Elements(1:columns_Number,3)=Nodes_cordinates1(Bays_Number+2:end,1);
Elements(1:columns_Number,4)=EAc*ones(columns_Number,1);
Elements(1:columns_Number,5)=EIc*ones(columns_Number,1);
%beams
BN=Nodes_cordinates1(Bays_Number+2:end,1);
BN2=BN;BN2(Bays_Number+1:Bays_Number+1:end)=[];
BN3=BN;BN3(1:Bays_Number+1:end)=[];
Elements(columns_Number+1:end-foundation_member,2)=BN2;
Elements(columns_Number+1:end-foundation_member,3)=BN3;
Elements(columns_Number+1:end-foundation_member,4)=EAb*ones(Beams_Number,1);
Elements(columns_Number+1:end-foundation_member,5)=EIb*ones(Beams_Number,1);
% foundation members
if foundn==1;
    Elements(columns_Number+Beams_Number+1:end,2)=supports(1:end-1,1);
    Elements(columns_Number+Beams_Number+1:end,3)=supports(2:end,1);
    Elements(columns_Number+Beams_Number+1:end,4)=EAf*ones(foundation_member,1);
    Elements(columns_Number+Beams_Number+1:end,5)=EIf*ones(foundation_member,1);
end
%    
Element_nodes=Elements(:,2:3);      % Elements connectivity
EA_mat=Elements(:,4);               % EA 
EI_mat=Elements(:,5);               % EI


if iip==1
    if iidof==1
% Checking Geometry
% % choice=framing_plot(Nodes_cordinates,Element_nodes);%if choice==1; break; end;
   
    end
end
%% Stiffness matrix calculation
Elements_number=length(Element_nodes)*1; %No discritization

%coordinates
xx=Nodes_cordinates(:,1);
yy=Nodes_cordinates(:,2);
Nodes_number=size(Nodes_cordinates,1);

% degrees of freedom
GDof= 3*Nodes_number;               %Global number of degrees of freedom

% Stiffness matrix & member lengths
[KG,Le]=stiffness2Dframe(GDof,Elements_number,Element_nodes,xx,yy,EI_mat,EA_mat);

%% Loading (only UDL on Beams)
% wb=-5000; %UDL on beams N/m

UDL=zeros(length(Element_nodes),1)';
%assign loading: [Beam No.1 UDL1; Beam No.2 UDL2; Beam No.3 UDL3; ....]
UDL_ass= [(columns_Number+1:Elements_Number);wb+zeros(1,length(columns_Number+1:Elements_Number))]'; % downword is -ve
UDL(UDL_ass(:,1))=UDL_ass(:,2);

% FEM
U=zeros(GDof,1);                    % Displacement vector
F=zeros(GDof,1);                    % External Force vector
FEMudl=zeros(GDof,1);               % Fixed end moment due to UDL
% UDL in members (counter-clockwise FEM is +ve)

for jj=1:length(UDL_ass(:,1))
    counter3=UDL_ass(jj,1);
    FEM_udl=zeros(GDof,1);
    FEM_udl(3*Element_nodes(counter3,1)-1)=-UDL(counter3)*Le(counter3)/2;
    FEM_udl(3*Element_nodes(counter3,2)-1)=-UDL(counter3)*Le(counter3)/2;
    FEM_udl(3*Element_nodes(counter3,1))=-UDL(counter3)*Le(counter3)^2/12;
    FEM_udl(3*Element_nodes(counter3,2))=UDL(counter3)*Le(counter3)^2/12;
    
    FEMudl=FEMudl+FEM_udl;
end
FEM=FEMudl;%+FEMF;

%% Supports settlement
supportDofs=sort([3*supports(:,1)-2; 3*supports(:,1)-1; 3*supports(:,1)]);
zeroDofs=find(supports(:,2:end)'==0); %Restrained Dofs

%assign settlment: 
settlement=reshape(support_settlement(:,2:end)',1,[])'; %vertical settlemnt at support(from LVDT) 

%% Solution
%solve for displacement F=KG*U+FEM
activeDof=setdiff(1:GDof,zeroDofs);
U(supportDofs)=settlement;
U(activeDof)=KG(activeDof,activeDof)\(F(activeDof)-FEM(activeDof)-KG(activeDof,supportDofs)*settlement);
Uxyz(:,2:4)=reshape(U,3,[])';Uxyz(:,1)=1:Nodes_Number; % Displacment in form [Node_No. ux, uy, uz, ...]

% Solve for Reactions > FL=KG*U+FEM;
FL=KG*U+FEM;                % Force in all nodes

Reactions=FL(zeroDofs);     % Reactions


end

