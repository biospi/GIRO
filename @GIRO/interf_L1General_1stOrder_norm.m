function [ctn, de_dCN] = interf_L1General_1stOrder_norm(OBJ_GIRO, CN)

stepSize = .001; 
%
% M_D = exp(OBJ_GIRO.BsplDictNorm * CN)'; 
%     
% % Log-Anscombe mean:
% de_dI = log((OBJ_GIRO.RefNorm - OBJ_GIRO.TarNorm .* M_D) + .375) + 2*log(2);
% 
% dI_dg = - OBJ_GIRO.TarNorm .* M_D ./ (OBJ_GIRO.RefNorm - OBJ_GIRO.TarNorm .* M_D + .375);
% 
% de_dCN = stepSize * ((de_dI .* dI_dg) * OBJ_GIRO.BsplDictNorm)' .* OBJ_GIRO.indCN;
% 
% ctn = norm((log((OBJ_GIRO.RefNorm - OBJ_GIRO.TarNorm .* M_D) + .375) + log(2)), 2)
 
M_D = exp(OBJ_GIRO.BsplDictNorm * CN)'; 

de_dI = OBJ_GIRO.RefNorm - log(Anscombe(OBJ_GIRO.TarNorm .* M_D));

dI_dg = - (OBJ_GIRO.TarNorm .* M_D) ./ (OBJ_GIRO.TarNorm .* M_D + .375);

de_dCN = stepSize * ((de_dI .* dI_dg) * OBJ_GIRO.BsplDictNorm)' .* OBJ_GIRO.indCN;

ctn = norm(OBJ_GIRO.RefNorm - log(Anscombe(OBJ_GIRO.TarNorm .* M_D)), 2) 




