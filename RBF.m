                                              
function [Rm]=RBF(x,fx,xi,dmj,delta)
N=length(x);

for i=1:N    
            rj=abs(x(i)-xi)/dmj;
            Rm(1,i)=((1-(rj/delta))^6)*(3+(18*rj/delta)+(35*(rj^2)/(delta^2)));     %RBF
                      
end

