function y=Ai(Xm,Ym,Zm,X,Y,Z)
y=[1 0 0 0        (Z-Zm)    -(Y-Ym)  ;
   0 1 0 -(Z-Zm)   0        (X-Xm)   ;
   0 0 1 (Y-Ym)    -(X-Xm)  0        ;
   0 0 0 1         0        0        ;
   0 0 0 0         1        0        ;
   0 0 0 0         0        1];       
