function [ ux,uz ]=u_KOR(z,x,So,deltaS,L,d,Vlt)  


uz= Vlt*(So+deltaS.*z/L);
ux=0.*x;


