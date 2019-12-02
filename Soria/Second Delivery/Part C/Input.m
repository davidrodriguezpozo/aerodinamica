function datos = Input(Vx,Vy)

    datos.Nx = Vx+2;
    datos.Ny = Vy+2;
    
    datos.Vx = Vx;
    datos.Vy = Vy;
    
    datos.L = 1;
    datos.H = 1;
    
    datos.uniform = true;
    datos.gamma = 1;
    
    datos.malla = 2;
    
    datos.dt_inicial = 0.001;
    datos.mu = 10e5;
    datos.time_final = 10;
    datos.delta = 10e-5;
    
    