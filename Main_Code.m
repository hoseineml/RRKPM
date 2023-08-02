clc;
clear;

x=linspace(-1,1,21);       %uniform arrange
fx=((25*x.^2)+1).^-1;


n=length(x);                      %number of node
np=(n-1)*10 ;                     %number of point


delta=4.1;                            % shape parameter for RBF
dmj=0.20;                              %radius of the supporting domain at nodes  For RBF 
h=1.15;                                %scaling parameter for RKPM

 for j=1:np
     
        xi=x(1)+(j-1)*((max(x)-min(x))/(np-1)); 

        % define Support domain for each point xi 
        xx=[];
        k=1;
        for i=1:n
            if  x(i)<xi+h && x(i)>xi-h %
                xx(k)=x(i);                                  
                k=k+1;  
            end
        end
         N=length(xx);
         fxx=((25*xx.^2)+1).^-1;
        
        % CALC R matrix
        Rm=zeros(N,N);
        Rn=zeros(N,n);
        
        
        for i=1:N
        Rm(i,:)=RBF(xx,fxx,xx(i),dmj,delta);
        Rn(i,:)=RKPM(x,fx,h,xx(i)); 
        RT=[Rm Rn];  
        end
        
        
        % CALC R(x) matrix
        Rmx=zeros(1,N);
        Rnx=zeros(1,n);
        Rmx=RBF(xx,fxx,xi,dmj,delta);
        [Rnx,W]=RKPMX(x,fx,h,xi,xx); 
        
        RTx=[Rmx Rnx];

        
        cond(RT);
        a=RT'*W*RT;
        cond(a);
        b=RT'*W;   
        cond(b);
        phix =RTx*(a\b) ;     %shape functions
        uhx(1,j)=phix*fxx' ;    
        uhx(2,j)=xi ; 
        error(j)=(uhx(1,j)-((25*xi.^2)+1).^-1);
        
 end

xxx=linspace(-1,1,100);
fxxx=((25*xxx.^2)+1).^-1;
figure(1)
for i=1:np
scatter(uhx(2,i),uhx(1,i),25);
scatter(x,fx,50,'filled','b');
plot(xxx,fxxx,'red','LineWidth',1.1);
hold on
end
legend('exact solution','RRKPM','known point','fontsize',13)
ylabel('uhx','fontsize',15); xlabel('x','fontsize',15)

figure(2)
for i=1:np
scatter(uhx(2,i),error(i),10);
hold on
end
ylabel('error','fontsize',15); xlabel('x','fontsize',15)


  