function error = ERROR(numeric,analytic)
    
    error = 0;

    for i = 1:size(numeric,1)
       for j = 1:size(numeric,2)
           %if analytic(i,j) < 0.0001 ||  analytic(i,j) > -0.0001
            %    a = abs(analytic(i,j)-numeric(i,j));
             %   disp('Algo va mal');
           %else 
               a = abs(analytic(i,j)-numeric(i,j));
           %end
            if a > error 
                error = a;
            end
       end
    end