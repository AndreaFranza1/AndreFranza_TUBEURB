function [Kelem] = KBern3D(E,d,dz,alfa,beta)

%==========================================================================
% Function per la generazione della matrice di rigidezza per l'elemento
% trave di Eulero-Bernoulli
%--------------------------------------------------------------------------
% Dati di input:
% Xi, Xf       	->  Coordinate intrinseche del nodo iniziale e del nodo 
%                   finale dell'elemento;
% E             ->  Modulo elastico del materiale costituente l'elemento
%                   U.M.    N/mm^2
% G             ->  Modulo tagliante del materiale costituente l'elemento
%                   U.M.    N/mm^2
% A             ->  Area della sezione trasversale dell'elemento
%                   U.M.    m^2
% I11           ->  Momento d'inerzia della sezione trasversale rispetto
%                   all'asse locale 1-1
%                   U.M.    m^4
% I22           ->  Momento d'inerzia della sezione trasversale rispetto
%                   all'asse locale 2-2
%                   U.M.    m^4
% I33           ->  Momento d'inerzia della sezione trasversale rispetto
%                   all'asse locale 3-3
%                   U.M.    m^4
%--------------------------------------------------------------------------
% Dati di output:
% 3DKBern       -> 	Matrice di rigidezza dell'elemento rispetto al sistema
%                   di riferimento globale
% Campo di spostamenti: (Piano 1-3)
%
%                             u3i                u3f
%                       (r22i |                  | (r22f
%                       u1i ->o------------------o-> u1f
%         3 |                 i                  f
%           |_ _
%               1
% Campo di spostamenti: (Piano 1-2)
%                             u2i                u2f
%                       )r33i |                  | )r33f
%                       u1i ->o------------------o-> u1f
%         2 |                 i                  f
%           |_ _
%               1
%
% Gli spostamenti nodali le rotazioni sono organizzate:
%
% U = [u1i, u2i, u3i, fi2i, fi3i,|u1f, u2f, u3f, fi2f, fi3f]
%
% Non si considera la componente rotazionale (torsionale) attorno l'asse 1
%--------------------------------------------------------------------------

% Assegnazione delle caratteristiche geometriche della sezione
A = pi.*d^2/4;
I11 = pi.*d^4/32;
I22 = pi.*d^4/64;
I33 = pi.*d^4/64;

G = E/(2*(1+0.2));

% Calcolo della lunghezza dell'elemento

Xi=[0 0 0];
Xf=[dz*tan(pi./180*beta) dz*tan(pi./180*alfa) dz];
L  = norm(Xf-Xi);

% Vettore di riferimento per la valutazione degli assi locali
% Non esiste una direzione preferenziale verso l'alto per un sistema di 
% coordinate locale. Tuttavia i sistemi di coordinate locali per nodi ed 
% elementi sono definiti rispetto alla direzione verso l'alto del sistema 
% globale, +Z. Potrebbero nascere dei problemi se l'asse dell'elemento
% coincide con l'asse Z; In questo caso i valori vengono inseriti
% manualmente.

Z = [0,0,1];

if abs(Xf-Xi)/L == abs(Z);
    % Caso in cui (Xf-Xi) coincide con l'asse Z
    x1 = (Xf-Xi)/norm(Xf-Xi);                              
    x3 = [0,1,0];           
    x2 = cross(x3,x1)/norm(cross(x3,x1));
else

    % Valutazione degli assi locali dell'elemento finito, in funzione della
    % convenzione adottata al punto precedente

    x1 = (Xf-Xi)/norm(Xf-Xi);
    x2 = cross(Z,x1)/norm(cross(Z,x1));
    x3 = cross(x1,x2)/norm(cross(x1,x2));
end

% Costruzione della Matrice delle trasformazioni degli assi (ROTAZIONI)
% Coseni direttori:
% lx1 = x1*[1;0;0];   lx2 = x2*[1;0;0];   lx3 = x3*[1;0;0];
% mx1 = x1*[0;1;0];   mx2 = x2*[0;1;0];   mx3 = x3*[0;1;0];
% nx1 = x1*[0;0;1];   nx2 = x2*[0;0;1];   nx3 = x3*[0;0;1];
% 
% Rot1 = [ lx1 , mx1 , nx1 ;
%          lx2 , mx2 , nx2 ;  
%          lx3 , mx3 , nx3 ]

Rot =  [ x1*[1;0;0], x1*[0;1;0], x1*[0;0;1], ;
         x2*[1;0;0], x2*[0;1;0], x2*[0;0;1], ;
         x3*[1;0;0], x3*[0;1;0], x3*[0;0;1], ];

MatRot = [       Rot, zeros(3,3), zeros(3,3), zeros(3,3);
          zeros(3,3),        Rot, zeros(3,3), zeros(3,3);
          zeros(3,3), zeros(3,3),        Rot, zeros(3,3);
          zeros(3,3), zeros(3,3), zeros(3,3),        Rot ];
        

% Matrice di rigidezza del singolo elemento, nel sistema di riferimento
% LOCALE x1,x2,x3

K = [ E*A/L,             0,             0,            0,            0,            0, -E*A/L,             0,             0,            0,             0,            0;
          0,  12*E*I33/L^3,             0,            0,            0,  6*E*I33/L^2,      0, -12*E*I33/L^3,             0,            0,             0,  6*E*I33/L^2;
          0,             0,  12*E*I22/L^3,            0, -6*E*I22/L^2,            0,      0,             0, -12*E*I22/L^3,            0,  -6*E*I22/L^2,            0;
          0,             0,             0,      G*I11/L/1e10,            0,            0,      0,             0,             0,     -G*I11/L/1e10,             0,            0;
          0,             0,  -6*E*I22/L^2,            0,    4*E*I22/L,            0,      0,             0,   6*E*I22/L^2,            0,     2*E*I22/L,            0;
          0,   6*E*I33/L^2,             0,            0,            0,    4*E*I33/L,      0,  -6*E*I33/L^2,             0,            0,             0,    2*E*I33/L;    
     -E*A/L,             0,             0,            0,            0,            0,  E*A/L,             0,             0,            0,             0,            0;
          0, -12*E*I33/L^3,             0,            0,            0, -6*E*I33/L^2,      0,  12*E*I33/L^3,             0,            0,             0, -6*E*I33/L^2;
          0,             0, -12*E*I22/L^3,            0,  6*E*I22/L^2,            0,      0,             0,  12*E*I22/L^3,            0,   6*E*I22/L^2,            0;
          0,             0,             0,     -G*I11/L/1e10,            0,            0,      0,             0,             0,      G*I11/L/1e10,             0,            0;      
          0,             0,  -6*E*I22/L^2,            0,    2*E*I22/L,            0,      0,             0,   6*E*I22/L^2,            0,     4*E*I22/L,            0;
          0,   6*E*I33/L^2,             0,            0,            0,    2*E*I33/L,      0,  -6*E*I33/L^2,             0,            0,             0,     4*E*I33/L ];    

     
% Matrice di rigidezza del singolo elemento, nel sistema di riferimento
% GLOBALE x,y,z

Kelem = MatRot.'*K*MatRot;

return

