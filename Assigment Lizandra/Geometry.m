function X = Geometry(div,type)

if strcmp(type,'wing') == 1
    switch div 
        case 16
            X = importdata ('NACA_2412_N_16.txt');
        case 32
            X = importdata ('NACA_2412_N_32.txt');    
        case 64
            X = importdata ('NACA_2412_N_64.txt');
        case 128
            X = importdata ('NACA_2412_N_128.txt');
        case 256
            X = importdata ('NACA_2412_N_256.txt');
        case 512
            X = importdata ('NACA_2412_N_512.txt');
    end
end

if strcmp(type,'tail') == 1
    switch div 
        case 16
            X = importdata ('NACA_0009_N_16.txt');
        case 32
            X = importdata ('NACA_0009_N_32.txt');    
        case 64
            X = importdata ('NACA_0009_N_64.txt');
        case 128
            X = importdata ('NACA_0009_N_128.txt');
        case 256
            X = importdata ('NACA_0009_N_256.txt');
        case 512
            X = importdata ('NACA_0009_N_512.txt');
    end
end