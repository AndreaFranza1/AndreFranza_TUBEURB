%% IMPORT INPUT DATA

% Settaggio dei fogli di lavoro, e dei percorsi
% Dati foglio di input

[blank, SYSname] = system('hostname');
if strcmp(SYSname(1:end-1),'PELKAYAM')
    
elseif strcmp(SYSname(1:end-1),'Andrea-UK')
folder = 'C:\Users\Andrea\Dropbox (Personal)\UoN\MATLAB\tunnel_pile_structure_int_finiteelement(TUST2017)\';
addpath('C:\Users\Andrea\Dropbox (Personal)\UoN\MATLAB\tunnel_pile_structure_int_finiteelement(TUST2017)\Script_Multilayerv2')
   

elseif strcmp(SYSname(1:end-1),'MSI')
folder = 'C:\Users\andre\Dropbox (Personal)\UoN\MATLAB\tunnel_pile_structure_int_finiteelement(TUST2017)\';
addpath('C:\Users\andre\Dropbox (Personal)\UoN\MATLAB\tunnel_pile_structure_int_finiteelement(TUST2017)\Script_Multilayerv2')

else
folder = 'C:\Users\Administrator\Dropbox (Personal)\UoN\MATLAB\tunnel_pile_structure_int_finiteelement(TUST2017)\';
addpath('C:\Users\Administrator\Dropbox (Personal)\UoN\MATLAB\tunnel_pile_structure_int_finiteelement(TUST2017)\Script_Multilayerv2')
   
end





ucap=xlsread(strcat(folder,namefileinp),sheetinp);

[ucap_pl,txt_pl,raw_pl]=xlsread(strcat(folder,namefileinp),sheetinp2);

%% READING INPUT DATA

%PROPERTIES OF THE SOIL
Es=ucap(:,1)';                  %[N/m2] Young's Modulus soil at pile shaft
Es=Es(~isnan(Es));
nis=ucap(:,2)';                 %[-] Poisson's ratio soil at pile shaft
nis=nis(~isnan(nis));
Eb=ucap(:,3)';                  %[N/m2] Young's Modulus soil at pile base
Eb=Eb(~isnan(Eb));
nib=ucap(:,4)';                 %[-] Poisson's ratio soil at pile base
nib=nib(~isnan(nib));


%GEOMETRY OF PILE GROUP
h_el=ucap(:,6)';                %h_el length finite element
h_el=h_el(~isnan(h_el));
X=ucap(:,7)';                   %X coordinates of pile heads
X=X(~isnan(X));
Y=ucap(:,8)';                   %Y coordinates of pile heads
Y=Y(~isnan(Y));
carat=ucap(:,12)';              %Vector Pile Properties
carat=carat(~isnan(carat));
alfa=ucap(:,9)';                %Inclination angle alfa of piles 
alfa=alfa(~isnan(alfa));
beta=ucap(:,10)';               %Inclination angle beta of piles 
beta=beta(~isnan(beta));

%PILE PROPERTIES
Ep=carat(1);                    %[N/m2] Young's Modulus Pile
d=carat(3);                     %[m]    Pile Diameter
L=carat(2);                     %[m]    Pile Length
Ip=pi*(d/2)^4/4;
Ap=pi*(d/2)^2;

%TUNNEL PROPERTIES
carat_tunnel=ucap(:,14)';             %Vector Tunnel Properties
carat_tunnel=carat_tunnel(~isnan(carat_tunnel));
ht=carat_tunnel(1);
Rt=carat_tunnel(2);
Vltp=carat_tunnel(3);
% epsilon=Vltp/100/2;
% epsilon=Vltp/100/2;
epsilon=(2-(4-4*Vltp/100).^0.5)./2;

%COORDINATE DEL NODO MASTER
Xm=carat(4);                    % cambia a seconda della palificata
Ym=carat(5);                    % cambia a seconda della palificata
Zm=carat(6);

% %Input Superstructure (beam or frame)
prop_frame=ucap(:,16)';              %vector frame properties
prop_frame=prop_frame(~isnan(prop_frame));
foundn=prop_frame(1); % Beam at z=0? yes:1;no:0;
Eframe=prop_frame(2);    %Pa  %30,600
num_storey=prop_frame(3);
spanx_frame=prop_frame(4);
spanz_frame=prop_frame(5);
bc=prop_frame(6); %m column
dc=prop_frame(7); %m column
bb=prop_frame(8); %m beam
db=prop_frame(9); %m beam
bf=prop_frame(10); %m foundation beam
df=prop_frame(11); %m foundation beam
npali=length(X);
num_bay=npali-1;

%%
% additionnal input % CHEN 
DH=ucap(:,5)';                 %[m] Thickness of each layer of soil 
DH=DH(~isnan(DH));             
DHS=DH;           

Ht=cumsum(DH);               %[m] The elevation of each layer of soil
Ht=[0,Ht(1:end-1)];

MSOILLS=length(DH);            % [-] Number of soil layers



%%



prop_pile  =ucap_pl(:,1);              
prop_pile  =prop_pile(~isnan(prop_pile));
SF0         =prop_pile(1);                  
Rs          =prop_pile(2);                  
Rb          =prop_pile(3);  



prop_alpha  =ucap_pl(:,4)';              
prop_alpha  =prop_alpha(~isnan(prop_alpha));
cu0         =prop_alpha(1); 
cuslope     =prop_alpha(2); 
alpha_cu    =prop_alpha(3); 
qb_lim_alpha      =prop_alpha(4); 

prop_beta   =ucap_pl(:,7)';              
prop_beta   =prop_beta(~isnan(prop_beta));
k0_coeff    =prop_beta(1); 
gammap      =prop_beta(2); 
beta_flim        =prop_beta(3); 
qb_lim_beta      =prop_beta(4); 
c_beta      =prop_beta(5); 

inputarray  =txt_pl(2:end,11)';  
 inputarray=inputarray(~cellfun('isempty',inputarray));