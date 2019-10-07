function Acarl = Poisson (datos, C)

Acarl = zeros(datos.Vx,datos.Vy);

for i = 2:datos.Nx-1
    for j = 2:datos.Ny-1
        
        coef_point_i = i;
        coef_point_j = j;
        coef_east_i = i+1;
        coef_east_j = j;
        coef_west_i = i-1;
        coef_west_j = j;
        coef_north_i = i;
        coef_north_j = j+1;
        coef_south_i = i;
        coef_south_j = j-1;
        
        coef_east(i-1,j-1) = IndexMatrix (datos.Nx,datos.Ny,coef_east_i,coef_east_j);
        coef_west(i-1,j-1) = IndexMatrix (datos.Nx,datos.Ny,coef_west_i,coef_west_j);
        coef_north(i-1,j-1) = IndexMatrix (datos.Nx,datos.Ny,coef_north_i,coef_north_j);
        coef_south(i-1,j-1) = IndexMatrix (datos.Nx,datos.Ny,coef_south_i,coef_south_j);
        coef_point(i-1,j-1) = IndexMatrix (datos.Nx,datos.Ny,coef_point_i,coef_point_j);
        
        delta_x = datos.L/datos.Nx;
        delta_y = datos.H/datos.Ny;
        
        if i==1 || j==1 || j == datos.Ny-1 || i == datos.Nx-1
             
        else
            
            coef_p = coef_point - datos.Nx;
            coef_e = coef_east - datos.Nx;
            coef_w = coef_west - datos.Nx;
            coef_n = coef_north - datos.Nx;
            coef_s = coef_south - datos.Nx;
        
%        if coef_west_i == 0
%            A (coef_west_i,coef_west_j) = 0;
%        else
            aw =  delta_x/delta_y;
            Acarl (coef_point(i-1,j-1), coef_west(i-1,j-1) ) =  aw; 
            Acarl (coef_west(i-1,j-1) , coef_point(i-1,j-1)) =  aw; 
            %Acarl2 ( coef_p , coef_w ) =  aw; 
%        end
        
%        if coef_east_i == datos.Nx
%            A (coef_east_i,coef_east_j) = 0;
%        else
            ae =  delta_x/delta_y;
            Acarl (coef_point(i-1,j-1), coef_east(i-1,j-1) ) = ae; 
            Acarl (coef_east(i-1,j-1) , coef_point(i-1,j-1) ) = ae; 
            %Acarl2 ( coef_p , coef_e ) =  ae; 
%        end
 
 
%        if coef_north_j == datos.Ny
%            A (coef_north_i,coef_north_j) = 0;
%        else
            an =  delta_x/delta_y;
            Acarl (coef_point(i-1,j-1), coef_north(i-1,j-1) ) = an; 
            Acarl ( coef_north(i-1,j-1) , coef_point(i-1,j-1) ) = an; 
            %Acarl2 ( coef_p , coef_n ) =  an; 
%        end
 
%        if coef_south_j == 0
%            A (coef_south_i,coef_south_j) = 0;
%        else
            as =  delta_x/delta_y;
            Acarl (coef_point(i-1,j-1), coef_south(i-1,j-1) ) = as; 
            Acarl (coef_south(i-1,j-1) , coef_point(i-1,j-1) ) = as; 
            %Acarl2 ( coef_p , coef_s ) =  as; 
%        end

             Acarl (coef_point(i-1,j-1),  coef_point(i-1,j-1) ) = -(ae + as + an + ae);
             %Acarl2 ( coef_p , coef_p ) =  ap; 
         end
    end
end

A_defi(1,:) = Acarl(7:9,2);
A_defi(2,:)  = Acarl(2:4,3);
A_defi(3,:)  = Acarl(2:4,4);

disp('Acabado');




function index = IndexMatrix (Nx,Ny,i,j)

    index = i + (j-1)*(Nx);