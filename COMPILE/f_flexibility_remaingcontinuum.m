function [FLEX_remaining,Xi_nod,Yi_nod,Zi_nod] = f_flexibility_remaingcontinuum(X,Y,z,Es,nis,nnodes,nelementi,d,X_additional)
%FLEXIBILITY MATRIX for the continuum





Y_additional=0*X_additional;
Y_all=[Y,Y_additional];
X_all=[X,X_additional];
npali_additional=length(X_all);

FLEX_remaining=zeros(npali_additional*(nelementi+1)*3);
nint=8*4;

for i=1:nnodes
	for ip=1:npali_additional
       z_global(i,ip)=z(i);
       x_global(i,ip)=X_all(ip);
       y_global(i,ip)=Y_all(ip);
	end
end


for i=1:npali_additional
    znod_tot(1+(nelementi+1)*(i-1):(nelementi+1)+(nelementi+1)*(i-1),1)=z';
    xnod_tot(1+(nelementi+1)*(i-1):1:(nelementi+1)+(nelementi+1)*(i-1),1)=X_all(i) ;
    ynod_tot(1+(nelementi+1)*(i-1):1:(nelementi+1)+(nelementi+1)*(i-1),1)=Y_all(i);
   
    Es_nod(1+(nelementi+1)*(i-1):(nelementi+1)+(nelementi+1)*(i-1),1)=Es(:,1);
    ni_nod(1+(nelementi+1)*(i-1):1:(nelementi+1)+(nelementi+1)*(i-1),1)=nis;
    

    
end


Gs_nod=Es_nod/2/(1+nis);




Zi_nod=znod_tot;
Xi_nod=xnod_tot;
Yi_nod=ynod_tot;

for jg=1:npali_additional*(nelementi+1)
    
 
    for ii_int=1:nint


        
%     Zj_nod(1:1:npali*(nelementi+1),1)=Zi_nod(jg);
%     Xj_nod(1:1:npali*(nelementi+1),1)=Xi_nod(jg)+d/2*cos(2*pi*(ii_int-1)/nint);
%     Yj_nod(1:1:npali*(nelementi+1),1)=Yi_nod(jg)+d/2*sin(2*pi*(ii_int-1)/nint);
% 

    
    if ii_int<=nint/2
        
    Zj_nod(1:1:npali_additional*(nelementi+1),1)=Zi_nod(jg);
    Xj_nod(1:1:npali_additional*(nelementi+1),1)=Xi_nod(jg)+d/2*cos(2*pi*(ii_int-1)/nint*2);
    Yj_nod(1:1:npali_additional*(nelementi+1),1)=Yi_nod(jg)+d/2*sin(2*pi*(ii_int-1)/nint*2);

    else
    Zj_nod(1:1:npali_additional*(nelementi+1),1)=Zi_nod(jg);
    Xj_nod(1:1:npali_additional*(nelementi+1),1)=Xi_nod(jg)+d/4*cos(2*pi*((ii_int-nint/2)-1)/nint*2);
    Yj_nod(1:1:npali_additional*(nelementi+1),1)=Yi_nod(jg)+d/4*sin(2*pi*((ii_int-nint/2)-1)/nint*2);
    
    end    
    

    
    
    [ uxij_1,uyij_1,uzij_1 ] = int_factor_mindlin_j_1_vect_CONT( Xi_nod,Yi_nod,Zi_nod,Xj_nod,Yj_nod,Zj_nod,Gs_nod,ni_nod);
    [ uxij_2,uyij_2,uzij_2 ] = int_factor_mindlin_j_2_vect_CONT( Xi_nod,Yi_nod,Zi_nod,Xj_nod,Yj_nod,Zj_nod,Gs_nod,ni_nod);
    [ uxij_3,uyij_3,uzij_3 ] = int_factor_mindlin_j_3_vect_CONT( Xi_nod,Yi_nod,Zi_nod,Xj_nod,Yj_nod,Zj_nod,Gs_nod,ni_nod);
    
    
    
    
    
    FLEX_remaining(1:3:npali_additional*(nelementi+1)*3,(jg-1)*3+1)=uxij_1/nint+FLEX_remaining(1:3:npali_additional*(nelementi+1)*3,(jg-1)*3+1);
    FLEX_remaining(2:3:npali_additional*(nelementi+1)*3,(jg-1)*3+1)=uyij_1/nint+FLEX_remaining(2:3:npali_additional*(nelementi+1)*3,(jg-1)*3+1);
    FLEX_remaining(3:3:npali_additional*(nelementi+1)*3,(jg-1)*3+1)=uzij_1/nint+FLEX_remaining(3:3:npali_additional*(nelementi+1)*3,(jg-1)*3+1);
    
    FLEX_remaining(1:3:npali_additional*(nelementi+1)*3,(jg-1)*3+2)=uxij_2/nint+FLEX_remaining(1:3:npali_additional*(nelementi+1)*3,(jg-1)*3+2);
    FLEX_remaining(2:3:npali_additional*(nelementi+1)*3,(jg-1)*3+2)=uyij_2/nint+FLEX_remaining(2:3:npali_additional*(nelementi+1)*3,(jg-1)*3+2);
    FLEX_remaining(3:3:npali_additional*(nelementi+1)*3,(jg-1)*3+2)=uzij_2/nint+FLEX_remaining(3:3:npali_additional*(nelementi+1)*3,(jg-1)*3+2);

    FLEX_remaining(1:3:npali_additional*(nelementi+1)*3,(jg-1)*3+3)=uxij_3/nint+FLEX_remaining(1:3:npali_additional*(nelementi+1)*3,(jg-1)*3+3);
    FLEX_remaining(2:3:npali_additional*(nelementi+1)*3,(jg-1)*3+3)=uyij_3/nint+FLEX_remaining(2:3:npali_additional*(nelementi+1)*3,(jg-1)*3+3);
    FLEX_remaining(3:3:npali_additional*(nelementi+1)*3,(jg-1)*3+3)=uzij_3/nint+FLEX_remaining(3:3:npali_additional*(nelementi+1)*3,(jg-1)*3+3);

    
    end
    
end


end

