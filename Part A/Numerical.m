function [conv_u diff_u conv_v diff_v] = Numerical (datos, C, u, v)

        conv_u = convection(u,v,datos.L,datos.H);
        diff_u = diffusion_u(u,datos.L);

        conv_v = convection(v,u,datos.H,datos.L);
        diff_v = diffusion_u(v,datos.H);
