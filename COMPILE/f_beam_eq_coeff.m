function [ EI_mat,EA_mat ] = f_beam_eq_coeff( H,Xheadp,Xheadp_i,Xheadp_f,Eframe,bf,df,bc,dc,bb,db,spanz_frame,spanx_frame,foundn,num_storey )
%F_BEAM_EQ_COEFF Summary of this function goes here
%   Detailed explanation goes here

num_p=length(Xheadp);

%I should neglect zone not in deflection outside tunnelling influence zone
if Xheadp_i<=(2.5*(-0.5*H))
    frXheadp_i=(2.5*(-0.5*H));
else 
    frXheadp_i=Xheadp_i;
end
if Xheadp_f>=(2.5*(-0.5*H))
    frXheadp_f=(2.5*(0.5*H));
else 
    frXheadp_f=Xheadp_i;
end

if frXheadp_i<=(-0.5*H)&&frXheadp_f<=(-0.5*H)
    num_bay_hogg_neg=(frXheadp_f-frXheadp_i)/spanx_frame;
    num_bay_sag=0;
    num_bay_hogg_pos=0;
elseif frXheadp_i<=(-0.5*H)&&frXheadp_f<=(0.5*H)&&frXheadp_f>=(-0.5*H)
    num_bay_hogg_neg=(-0.5*H-frXheadp_i)/spanx_frame;
    num_bay_sag=(frXheadp_f-0.5*H)/spanx_frame;
    num_bay_hogg_pos=0;
elseif frXheadp_i<=(-0.5*H)&&frXheadp_f>=(0.5*H)
    num_bay_hogg_neg=(-0.5*H-frXheadp_i)/spanx_frame;
    num_bay_sag=H/spanx_frame;
    num_bay_hogg_pos=(frXheadp_f-0.5*H)/spanx_frame;
elseif frXheadp_i>=(-0.5*H)&&frXheadp_i<=(0.5*H)&&frXheadp_f<=(0.5*H)
    num_bay_hogg_neg=0;
    num_bay_sag=(frXheadp_f-frXheadp_i)/spanx_frame;
    num_bay_hogg_pos=0;
elseif frXheadp_i>=(-0.5*H)&&frXheadp_i<=(0.5*H)&&frXheadp_f>=(0.5*H)
    num_bay_hogg_neg=0;
    num_bay_sag=(0.5*H-frXheadp_i)/spanx_frame;
    num_bay_hogg_pos=(frXheadp_f-0.5*H)/spanx_frame;
elseif frXheadp_i>=(0.5*H)
    num_bay_hogg_neg=0;
    num_bay_sag=0;
    num_bay_hogg_pos=(frXheadp_f-frXheadp_i)/spanx_frame;
end

klcfr=Eframe*(bc*dc^3/12)/spanz_frame;
kbfr=Eframe*(bb*db^3/12)/spanx_frame;

c_eq_fr_hogg_neg=1+num_bay_hogg_neg^2*(2*klcfr)/(2*klcfr+kbfr);
c_eq_fr_sag=1+num_bay_sag^2*(2*klcfr)/(2*klcfr+kbfr);
c_eq_fr_hogg_pos=1+num_bay_hogg_pos^2*(2*klcfr)/(2*klcfr+kbfr);


% %Equivalent square beam at foundation level modelling frame
if foundn==1
	EI_eq_fr_hogg_neg= Eframe*(bf*df^3/12+num_storey*c_eq_fr_hogg_neg*bb*db^3/12);
    EI_eq_fr_hogg_pos= Eframe*(bf*df^3/12+num_storey*c_eq_fr_hogg_pos*bb*db^3/12);
    EI_eq_fr_sag= Eframe*(bf*df^3/12+num_storey*c_eq_fr_sag*bb*db^3/12);
    EA_eq_fr_hogg_neg=Eframe*(bf*df+num_storey*bb*db);
    EA_eq_fr_hogg_pos=Eframe*(bf*df+num_storey*bb*db);
    EA_eq_fr_sag=Eframe*(bf*df+num_storey*bb*db);
else
	EI_eq_fr_hogg_neg= Eframe*(bf*df+num_storey*bb*db);(num_storey*c_eq_fr_hogg_neg*bb*db^3/12);
    EI_eq_fr_hogg_pos= Eframe*(bf*df+num_storey*bb*db);(num_storey*c_eq_fr_hogg_pos*bb*db^3/12);
    EI_eq_fr_sag= Eframe*(bf*df+num_storey*bb*db);(num_storey*c_eq_fr_sag*bb*db^3/12);
    EA_eq_fr_hogg_neg=Eframe*(num_storey*bb*db);
    EA_eq_fr_hogg_pos= Eframe*(num_storey*bb*db);
    EA_eq_fr_sag=Eframe*(num_storey*bb*db);
end

EI_mat=zeros(1,num_p);
EA_mat=zeros(1,num_p);


for ip=1:num_p
    if Xheadp(ip)>=(-0.5*H)&&Xheadp(ip)<=(0.5*H)
        EI_mat(ip)=EI_eq_fr_sag;
        EA_mat(ip)=EA_eq_fr_sag;
    elseif Xheadp(ip)<(-0.5*H)
        EI_mat(ip)=EI_eq_fr_hogg_neg;
        EA_mat(ip)=EA_eq_fr_hogg_neg;
    else
        EI_mat(ip)=EI_eq_fr_hogg_pos;
        EA_mat(ip)=EA_eq_fr_hogg_pos;
    end
end


end

