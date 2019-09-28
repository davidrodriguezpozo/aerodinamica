function main
    syms x y u v
    % just an example:
    u=sin(x)*cos(y)
    v=cos(x)*sin(y)
    diff(u,x) % du/dx
    diff(v,y) % dv/dy
    diff(u,x,x) % d^2u/dx^2
    field = diff(u,x)+diff(v,y)
    % now, if we want to evaluate numerically
    % the best is transform the symbolic expression
    % ’field' to a function
    class(field)
    f=matlabFunction(field, 'Vars', [x y]);
    % (the argument 'Vars' forces the function f +
    % to have two arguments, x and y, precisely in this order
    class(f)
    f(0.4,0.2)

    
function cu = convection (u,v,L,H)
    
    % CDS : Central Difference Scheme
        ue = (u(i+1,j)+u(i,j))/2;
        uw = (u(i+1,j)+u(i,j))/2;
        vn = (v(i,j+1)+u(i,j))/2;
        vs = (v(i,j-1)+u(i,j))/2;
        un = (u(i,j+1)+u(i,j))/2;
        us = (u(i,j-1)+u(i,j))/2;
   
   % Probar Upwind més endavant
        
        Nx = length(u,1);
        Ny = length(u,2);
        
        dx = L/N; dy = H/Ny;
   
   % Mass Fluxes 
   
        Fe = ue * dy;
        Fw = uw * dy;
        Fn = vn * dx;
        Fs = vs * dx;
    
   % Convection term
   
   cu = Fe*ue - Fw*uw + Fn*un - Fs*us;

   
    function [diff_x_anal diff_y_anal conv_x_anal conv_y_anal] = ANALYTIC (X, Y)
       %syms x y u v
       % just an example:
       
       
        u = datos.F*cos(2*pi*x)*sin(2*pi*y)
        v = -datos.F*cos(2*pi*y)*sin(2*pi*x)
        
        % du/dx = ?2*pi*sin(2*pi*y)*sin(2*pi*x)
        % dv/dy = 2*pi*sin(2*pi*x)*sin(2*pi*y)
        % d^2u/dx^2 = ?4*pi^2*sin(2*pi*y)*cos(2*pi*x)
        % d^2v/dy^2 = 4*pi^2*sin(2*pi*x)*cos(2*pi*y)
        
        % u * du/dx = datos.F^2 * ( cos(2*pi*x)*sin(2*pi*y) ) * ( ?2*pi*sin(2*pi*y)*sin(2*pi*x) ) 
        % v * du/dx = datos.F^2 * ( -cos(2*pi*y)*sin(2*pi*x) ) * ( ?2*pi*sin(2*pi*y)*sin(2*pi*x) )
        % u * dv/dx = datos.F^2 * ( cos(2*pi*x)*sin(2*pi*y) ) * ( 2*pi*sin(2*pi*x)*sin(2*pi*y) )
        % v * dv/dx = datos.F^2 * ( -cos(2*pi*y)*sin(2*pi*x) ) * ( 2*pi*sin(2*pi*x)*sin(2*pi*y) )
        
        diff(u,x) % du/dx
        diff(v,y) % dv/dy
        diff(u,x,x) % d^2u/dx^2
        
        exp_diff_x_anal = diff(u,x,x);
        exp_diff_y_anal = diff(v,y,y);
        exp_conv_x_anal = u*diff(u,x) + u*diff(v,x) ;
        exp_conv_y_anal = v*diff(u,x) + v*diff(v,x) ;
        
        Nx = size(X,1);
        Ny = size(X,2);
        
        % Calculem una matriu amb els valors analitics en cada punt de la
        % malla (malla del resultat a comprovar)
        
        for i=1:Nx
            for j=1:Ny
                
                xi = X(i,j);
                yi = Y(i,j);
                
                %   DIFFUSSIO : x
                    diff_x_anal(i,j) = sym2numeric(exp_diff_x_anal(xi,yi)); 
        
                %   DIFFUSSIO : y
                    diff_y_anal(i,j) = sym2numeric(exp_diff_y_anal(xi,yi)); 
        
                %   CONVECCIO : x
                    conv_x_anal(i,j) = sym2numeric(exp_diff_y_anal(xi,yi)); 
                    
                %   CONVECCIO : y
                    conv_y_anal(i,j) = sym2numeric(exp_conv_y_anal(xi,yi)); 
            
            end
        end             
             
   
   function error = ERROR(numeric,anaylitic)
       
       A = abs(analytic-numeric)
       
       error = max(A);
        
        
        
        