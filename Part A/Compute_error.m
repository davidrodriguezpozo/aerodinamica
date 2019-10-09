function [error_du error_cu error_dv error_cv] = Compute_error ( S, k) 

error_du = error(S.diffu,S.du_anal);
error_cu = error(S.cu,S.cu_anal);
error_dv = error(S.diffv,S.dv_anal);
error_cv = error(S.cv,S.cv_anal);

%disp('Error calculat');