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
error_du = error(S.diffu,S.du_anal);
error_cu = error(S.cu,S.cu_anal);
error_dv = error(S.diffv,S.dv_anal);
error_cv = error(S.cv,S.cv_anal);
