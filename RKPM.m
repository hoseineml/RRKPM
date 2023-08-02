function [Rnx]= RKPM ( x , fx , h , xi)

    n = length(x) ;    
    u = 0 ;
    Rnx=zeros(1,n);
    phi = zeros(n,1);
    W =   zeros(n,1) ;
    m0 = 0;
    m1 = 0;
    m2 = 0;
    
    %dx=zeros(1,n);
    
    for i=1:n
        if i<n
        dx(i)=x(i+1)-x(i);
        elseif i==n
            dx(i)=x(i)-x(i-1);
        end
    end   
            

    for i=1:n
        d = abs(xi - x(i))/h ;    

        if d<=1 
              %W(i)= (2/3)-(4*(d)^2)+(4*(d^3));                %cubic spline
              %W(i)= (1-(6*(d^2))+(8*(d^3))-(3*(d^4))) ;       %quartic spline
              W(i)= exp(-(d^2)/(0.3^2)) ;                     % exponential
               %nn=1.1;                                         % exponential gauss
               %W(i)= (1/((pi^(nn/2))*(h^nn)))*exp(-(d^2));     % exponential gauss
        
        else W(i) = 0 ;
        end
             
        m0 = m0 + W(i)*dx(i) ;
        m1 = m1 + ((xi-x(i)))*W(i)*dx(i) ;
        m2 = m2 + ((xi-x(i)))^2*W(i)*dx(i) ;
    end
    
    
        C1 = (m2)/(m0*m2-m1^2) ;
        C2 = -(m1)/(m0*m2-m1^2) ;
    
    for i=1:n
        Ch = C1 + C2*(xi-x(i)) ;
        %u = u + fx(i) * Ch * W(i) * dx ;
        phi(i) = Ch * W(i) * dx(i)  ;
        Rnx(i) = fx(i) * Ch * W(i) * dx(i) ;
    end

end

