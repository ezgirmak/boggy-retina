function C = CellComp4Layer_Coefficients(Coeff,eta_V,eta_F,eta_L,eta_I,sigma_V,xi_T_f,sigma_L_y,sigma_I,Yii, h_F, m)

if strcmpi(Coeff,'A1')
    
    C = 0;
    
elseif strcmpi(Coeff,'B1')
    
    C = (m.*exp(Yii.*eta_V + 2.*eta_F.*h_F).*(eta_V.*sigma_V - ...
    eta_F.*xi_T_f).*(eta_L.*sigma_L_y + eta_F.*xi_T_f) - ...
    m.*exp(Yii.*eta_V).*(eta_V.*sigma_V + eta_F.*xi_T_f).*(eta_L.*sigma_L_y - ...
    eta_F.*xi_T_f))./(eta_V);
    
elseif strcmp(Coeff,'B2')
    
    C = (eta_V.*m.*sigma_V.*cosh(Yii.*eta_V).*(eta_F.*xi_T_f - ...
    eta_L.*sigma_L_y + eta_L.*sigma_L_y.*exp(2.*eta_F.*h_F) + ...
    eta_F.*xi_T_f.*exp(2.*eta_F.*h_F)) - eta_F.^2.*m.*xi_T_f.^2.*sinh(Yii.*eta_V) ...
    + eta_F.*eta_L.*m.*sigma_L_y.*xi_T_f.*sinh(Yii.*eta_V) + ...
    eta_F.*m.*xi_T_f.*exp(2.*eta_F.*h_F).*sinh(Yii.*eta_V).*(eta_L.*sigma_L_y + ...
    eta_F.*xi_T_f))./(eta_V);
    
elseif strcmpi(Coeff,'C1')
    
    C = -(2.*m.*sigma_V.*exp(Yii.*eta_V).*(eta_L.*sigma_L_y - ...
    eta_F.*xi_T_f));
    
elseif strcmpi(Coeff,'C2')
    
    C = (2.*m.*sigma_V.*exp(Yii.*eta_V).*exp(2.*eta_F.*h_F).*(eta_L.*sigma_L_y + ...
    eta_F.*xi_T_f));
    
elseif strcmpi(Coeff,'D2')
    
    C =(4.*eta_F.*m.*sigma_V.*xi_T_f.*exp(h_F.*(eta_F + ...
    eta_L)).*exp(Yii.*eta_V));
    
end