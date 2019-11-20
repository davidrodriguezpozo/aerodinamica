function [R_u R_v R_u_prev R_v_prev] = ComputeR (datos, C, u, v)

[S.cu S.diffu S.cv S.diffv] = Numerical (datos, C, u, v);
