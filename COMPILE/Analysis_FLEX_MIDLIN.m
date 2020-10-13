clear all


tic

d=0.8;
h_el=1;
Lp=25;
npali=4;
nnodes=Lp/h_el+1;
nelementi=Lp/h_el;
Zi_nod=[	0	1	2	3	4	5	6	7	8	9	10	11	12	13	14	15	16	17	18	19	20	21	22	23	24	25	0	1	2	3	4	5	6	7	8	9	10	11	12	13	14	15	16	17	18	19	20	21	22	23	24	25	0	1	2	3	4	5	6	7	8	9	10	11	12	13	14	15	16	17	18	19	20	21	22	23	24	25	0	1	2	3	4	5	6	7	8	9	10	11	12	13	14	15	16	17	18	19	20	21	22	23	24	25]';
Xi_nod=	[4.5	4.5	4.5	4.5	4.5	4.5	4.5	4.5	4.5	4.5	4.5	4.5	4.5	4.5	4.5	4.5	4.5	4.5	4.5	4.5	4.5	4.5	4.5	4.5	4.5	4.5	6.9	6.9	6.9	6.9	6.9	6.9	6.9	6.9	6.9	6.9	6.9	6.9	6.9	6.9	6.9	6.9	6.9	6.9	6.9	6.9	6.9	6.9	6.9	6.9	6.9	6.9	4.5	4.5	4.5	4.5	4.5	4.5	4.5	4.5	4.5	4.5	4.5	4.5	4.5	4.5	4.5	4.5	4.5	4.5	4.5	4.5	4.5	4.5	4.5	4.5	4.5	4.5	6.9	6.9	6.9	6.9	6.9	6.9	6.9	6.9	6.9	6.9	6.9	6.9	6.9	6.9	6.9	6.9	6.9	6.9	6.9	6.9	6.9	6.9	6.9	6.9	6.9	6.9]';
Yi_nod=	[-1.2	-1.2	-1.2	-1.2	-1.2	-1.2	-1.2	-1.2	-1.2	-1.2	-1.2	-1.2	-1.2	-1.2	-1.2	-1.2	-1.2	-1.2	-1.2	-1.2	-1.2	-1.2	-1.2	-1.2	-1.2	-1.2	-1.2	-1.2	-1.2	-1.2	-1.2	-1.2	-1.2	-1.2	-1.2	-1.2	-1.2	-1.2	-1.2	-1.2	-1.2	-1.2	-1.2	-1.2	-1.2	-1.2	-1.2	-1.2	-1.2	-1.2	-1.2	-1.2	1.2	1.2	1.2	1.2	1.2	1.2	1.2	1.2	1.2	1.2	1.2	1.2	1.2	1.2	1.2	1.2	1.2	1.2	1.2	1.2	1.2	1.2	1.2	1.2	1.2	1.2	1.2	1.2	1.2	1.2	1.2	1.2	1.2	1.2	1.2	1.2	1.2	1.2	1.2	1.2	1.2	1.2	1.2	1.2	1.2	1.2	1.2	1.2	1.2	1.2	1.2	1.2]';
Gs_nod=	[8000000	8000000	8000000	8000000	8000000	8000000	8000000	8000000	8000000	8000000	8000000	8000000	8000000	8000000	8000000	8000000	8000000	8000000	8000000	8000000	8000000	8000000	8000000	8000000	8000000	8000000	8000000	8000000	8000000	8000000	8000000	8000000	8000000	8000000	8000000	8000000	8000000	8000000	8000000	8000000	8000000	8000000	8000000	8000000	8000000	8000000	8000000	8000000	8000000	8000000	8000000	8000000	8000000	8000000	8000000	8000000	8000000	8000000	8000000	8000000	8000000	8000000	8000000	8000000	8000000	8000000	8000000	8000000	8000000	8000000	8000000	8000000	8000000	8000000	8000000	8000000	8000000	8000000	8000000	8000000	8000000	8000000	8000000	8000000	8000000	8000000	8000000	8000000	8000000	8000000	8000000	8000000	8000000	8000000	8000000	8000000	8000000	8000000	8000000	8000000	8000000	8000000	8000000	8000000]';
ni_nod=	[0.5	0.5	0.5	0.5	0.5	0.5	0.5	0.5	0.5	0.5	0.5	0.5	0.5	0.5	0.5	0.5	0.5	0.5	0.5	0.5	0.5	0.5	0.5	0.5	0.5	0.5	0.5	0.5	0.5	0.5	0.5	0.5	0.5	0.5	0.5	0.5	0.5	0.5	0.5	0.5	0.5	0.5	0.5	0.5	0.5	0.5	0.5	0.5	0.5	0.5	0.5	0.5	0.5	0.5	0.5	0.5	0.5	0.5	0.5	0.5	0.5	0.5	0.5	0.5	0.5	0.5	0.5	0.5	0.5	0.5	0.5	0.5	0.5	0.5	0.5	0.5	0.5	0.5	0.5	0.5	0.5	0.5	0.5	0.5	0.5	0.5	0.5	0.5	0.5	0.5	0.5	0.5	0.5	0.5	0.5	0.5	0.5	0.5	0.5	0.5	0.5	0.5	0.5	0.5]';


%% Flexibility matrix for the continuum method



KKsoil=zeros(npali*(nelementi+1)*6);


FLEX_pipe=zeros(npali*(nelementi+1)*3);
nint=8*4;

%FLEX_pipeIBILITY MATRIX for the continuum
for jg=1:npali*(nelementi+1)
    
    for ii_int=1:nint

% 
%         
%     Zj_nod(1:1:npali*(nelementi+1),1)=Zi_nod(jg);
%     Xj_nod(1:1:npali*(nelementi+1),1)=Xi_nod(jg)+d/2*cos(2*pi*(ii_int-1)/nint);
%     Yj_nod(1:1:npali*(nelementi+1),1)=Yi_nod(jg)+d/2*sin(2*pi*(ii_int-1)/nint);


    
    if ii_int<=nint/2
        
    Zj_nod(1:1:npali*(nelementi+1),1)=Zi_nod(jg);
    Xj_nod(1:1:npali*(nelementi+1),1)=Xi_nod(jg)+d/2*cos(2*pi*(ii_int-1)/nint*2);
    Yj_nod(1:1:npali*(nelementi+1),1)=Yi_nod(jg)+d/2*sin(2*pi*(ii_int-1)/nint*2);

    else
    Zj_nod(1:1:npali*(nelementi+1),1)=Zi_nod(jg);
    Xj_nod(1:1:npali*(nelementi+1),1)=Xi_nod(jg)+d/4*cos(2*pi*((ii_int-nint/2)-1)/nint*2);
    Yj_nod(1:1:npali*(nelementi+1),1)=Yi_nod(jg)+d/4*sin(2*pi*((ii_int-nint/2)-1)/nint*2);
    
    end    
    
    dZj_nod(:,jg)=Zi_nod(:)-Zj_nod(:);
    dXj_nod(:,jg)=Xi_nod(:)-Xj_nod(:);
    dYj_nod(:,jg)=Yi_nod(:)-Yj_nod(:); 
    

    
    
        
    
    
    [ uxij_1,uyij_1,uzij_1 ] = int_factor_mindlin_j_1_vect_CONT( Xi_nod,Yi_nod,Zi_nod,Xj_nod,Yj_nod,Zj_nod,Gs_nod,ni_nod);
    [ uxij_2,uyij_2,uzij_2 ] = int_factor_mindlin_j_2_vect_CONT( Xi_nod,Yi_nod,Zi_nod,Xj_nod,Yj_nod,Zj_nod,Gs_nod,ni_nod);
    %I need to check the following function
    [ uxij_3,uyij_3,uzij_3 ] = int_factor_mindlin_j_3_vect_CONT( Xi_nod,Yi_nod,Zi_nod,Xj_nod,Yj_nod,Zj_nod,Gs_nod,ni_nod);
    
    
    
    
    
    FLEX_pipe(1:3:npali*(nelementi+1)*3,(jg-1)*3+1)=uxij_1/nint+FLEX_pipe(1:3:npali*(nelementi+1)*3,(jg-1)*3+1);
    FLEX_pipe(2:3:npali*(nelementi+1)*3,(jg-1)*3+1)=uyij_1/nint+FLEX_pipe(2:3:npali*(nelementi+1)*3,(jg-1)*3+1);
    FLEX_pipe(3:3:npali*(nelementi+1)*3,(jg-1)*3+1)=uzij_1/nint+FLEX_pipe(3:3:npali*(nelementi+1)*3,(jg-1)*3+1);
    
    FLEX_pipe(1:3:npali*(nelementi+1)*3,(jg-1)*3+2)=uxij_2/nint+FLEX_pipe(1:3:npali*(nelementi+1)*3,(jg-1)*3+2);
    FLEX_pipe(2:3:npali*(nelementi+1)*3,(jg-1)*3+2)=uyij_2/nint+FLEX_pipe(2:3:npali*(nelementi+1)*3,(jg-1)*3+2);
    FLEX_pipe(3:3:npali*(nelementi+1)*3,(jg-1)*3+2)=uzij_2/nint+FLEX_pipe(3:3:npali*(nelementi+1)*3,(jg-1)*3+2);

    FLEX_pipe(1:3:npali*(nelementi+1)*3,(jg-1)*3+3)=uxij_3/nint+FLEX_pipe(1:3:npali*(nelementi+1)*3,(jg-1)*3+3);
    FLEX_pipe(2:3:npali*(nelementi+1)*3,(jg-1)*3+3)=uyij_3/nint+FLEX_pipe(2:3:npali*(nelementi+1)*3,(jg-1)*3+3);
    FLEX_pipe(3:3:npali*(nelementi+1)*3,(jg-1)*3+3)=uzij_3/nint+FLEX_pipe(3:3:npali*(nelementi+1)*3,(jg-1)*3+3);

    
    
    end



    
end
toc