% This function to calculate the stiffness of 2D frame
function [KG,Le]= ...
    stiffness2Dframe(GDof,Elements_number,Element_nodes,xx,yy,EI_mat,EA_mat);

%Global stiffnes matrix
KG=zeros(GDof);
for counter2=1:Elements_number
    EA=EA_mat(counter2);
    EI=EI_mat(counter2);
    
    %Calculation of element length & angle of inclination
    x1=xx(Element_nodes(counter2,1));
    y1=yy(Element_nodes(counter2,1));
    x2=xx(Element_nodes(counter2,2));
    y2=yy(Element_nodes(counter2,2));
    
    le=sqrt((x2-x1)^2+(y2-y1)^2); %Length of the memeber
    Le(counter2)=le;
    C=(x2-x1)/le;                 % Cosine angle of inclination
    S=(y2-y1)/le;                 % Sine angle of inclination
    
    %Element stiffness
    Ke=[EA/le  ,0          ,0          ,-EA/le ,0          ,0;
        0      ,12*EI/le^3 ,6*EI/le^2  ,0      ,-12*EI/le^3 ,6*EI/le^2;
        0      ,6*EI/le^2  ,4*EI/le    ,0      ,-6*EI/le^2 ,2*EI/le;
        -EA/le ,0          ,0          ,EA/le  ,0          ,0;
        0      ,-12*EI/le^3 ,-6*EI/le^2 ,0      ,12*EI/le^3 ,-6*EI/le^2;
        0      ,6*EI/le^2  ,2*EI/le    ,0      ,-6*EI/le^2 ,4*EI/le];
    
    %Transformation matrix, 
        T = [C ,S ,0 ,0  ,0 ,0;...
            -S ,C ,0 ,0  ,0 ,0;...
            0  ,0 ,1 ,0  ,0 ,0;...
            0  ,0 ,0 ,C  ,S ,0;...
            0  ,0 ,0 ,-S ,C ,0;...
            0  ,0 ,0 ,0  ,0 ,1];
   
    %Selection of element's Dofs
    ENs=Element_nodes(counter2,:);  %Nodes of this element
    EDofs=[[1 2 3]+3*(ENs(1)-1) [1 2 3]+3*(ENs(2)-1)]; %Degress of freedom associated to this element
        
    %Assembly of KG
    KG(EDofs,EDofs)=KG(EDofs,EDofs)+T'*Ke*T;end
        
        

