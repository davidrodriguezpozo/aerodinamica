%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%       FUNCTION COMPUTE_ERROR      %%%%%%%%%%%%%
%%% This function computes, as it is said in the name,  %%%
%%% the error between the numeric and analitic solution %%%
%%% with the final purpose to validate the code of the  %%%
%%% Part A.                                             %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [error_du error_cu error_dv error_cv] = Compute_error ( S, k) 
% The Outputs: diffusive u term error, diffusive v term error, convective u term error, diffusive v term error.
% The Inputs: the solution (S) 
error_du = error_ind(S.diffu,S.du_analytic);
error_cu = error_ind(S.cu,S.cu_analytic);
error_dv = error_ind(S.diffv,S.dv_analytic);
error_cv = error_ind(S.cv,S.cv_analytic);
