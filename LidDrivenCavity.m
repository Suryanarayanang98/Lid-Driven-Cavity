%% Code for Lid Driven Cavity

% Geometry Set up
gridlength=0.05; % delta x = delta y = gridlength= 0.05
length=2; 
height=2;
Re=100; % Reynold's Number
U=3; % Lid Velocity
mu=U*length/Re;  % Coefficient of viscosity

% Carrying out Stability condition
a=gridlength*gridlength/(2*mu);
b=gridlength/U;
timestep=min(a,b);

%Computing Number of grids in x and y directions
ni=height/gridlength;   % 40
nj=length/gridlength;    % 40

%Initialising stream function and vorticity values
psi=zeros(ni,nj);
zeta=zeros(ni,nj);
zetaNext=zeros(ni,nj); % Vorticity values for the next timestep

%Residual Calculations
niter=zeros(1,12000);
errorpsimat=zeros(1,12000);
errorpsi=0;
errorzeta=0;
errorp=0;
errorzetamat=zeros(1,12000);

%Initialising Pressure ,Horizontal and vertical velocities
u=zeros(ni,nj);
v=zeros(ni,nj);
p=zeros(ni,nj);

% Initialisation of horizontal velocity at the lid
for i=1:ni-1
    u(1,i)=U;
end

for iter=1:120000
    %% Assigning zero to errors for each iteration
    errorpsi=0;
    errorzeta=0;
    %% Vorticity Transport Equation
    % Boundary Conditions
    for k=1:ni
        zeta(1,k)=2*(-gridlength*U -psi(2,k))/(gridlength*gridlength);  %Top Lid
        zeta(nj,k)=2*-psi(nj-1,k)/(gridlength*gridlength); %Bottom wall
    end
    for l=2:nj-1
        zeta(l,1)=2*-psi(l,2)/(gridlength*gridlength); % Left wall
        zeta(l,ni)=2*-psi(l,ni-1)/(gridlength*gridlength); % Right wall
    end     

    for i=2:(nj-1)  % ROW
        for j=2:(ni-1) %COLUMN
            RHS=((zeta(i,j+1)-2*zeta(i,j)+zeta(i,j-1))/(gridlength*gridlength) + (zeta(i-1,j)-2*zeta(i,j)+zeta(i+1,j))/(gridlength*gridlength));
            convective_terms= (u(i,j+1)*zeta(i,j+1)-u(i,j-1)*zeta(i,j-1))/(2*gridlength) + (v(i-1,j)*zeta(i-1,j)-v(i+1,j)*zeta(i+1,j))/(2*gridlength);
            zetaNext(i,j)=zeta(i,j)+timestep*( (RHS/Re) -convective_terms);
            errorzeta=errorzeta+abs(zetaNext(i,j)-zeta(i,j)); %Residual 
        end
    end

    %% Stream Function Equation
    for i=2:(nj-1)  % ROW
        for j=2:(ni-1) %COLUMN
            old=psi(i,j);
            psi(i,j)=(psi(i+1,j)+psi(i-1,j)+psi(i,j+1)+psi(i,j-1)+(gridlength*gridlength)*zetaNext(i,j))/4;
            errorpsi=errorpsi+abs(psi(i,j)-old); % Residual Calculation
        end
    end
    %Boundary condition updation
    for k=1:ni
        zetaNext(1,k)=2*(-gridlength*U -psi(2,k))/(gridlength*gridlength); % Top lid
        errorzeta=errorzeta+abs(zetaNext(1,k)-zeta(1,k)); % Residual
        zetaNext(nj,k)=2*-psi(nj-1,k)/(gridlength*gridlength); % Bottom wall
        errorzeta=errorzeta+abs(zetaNext(nj,k)-zeta(nj,k)); % Residual
    end
    for l=2:nj-1
        zetaNext(l,1)=2*-psi(l,2)/(gridlength*gridlength); % Left wall
        errorzeta=errorzeta+abs(zetaNext(l,1)-zeta(l,1)); % Residual
        zetaNext(l,ni)=2*-psi(l,ni-1)/(gridlength*gridlength); % Right wall
        errorzeta=errorzeta+abs(zetaNext(l,ni)-zeta(l,ni)); % Residual
    end     

    %% Getting Velocities
    for i=2:(nj-1)  % ROW
        for j=2:(ni-1) %COLUMN
            u(i,j)= (psi(i,j)-psi(i+1,j))/(1*gridlength); % Horizontal
            v(i,j)=-(psi(i,j)-psi(i,j-1))/(1*gridlength); %Vertical
        end
    end

    zeta=zetaNext; % Assigning zeta of next time step
    errorpsimat(1,iter)=errorpsi/(ni*nj); % Error in psi at each iteration
    errorzetamat(1,iter)=errorzeta/(ni*nj); % Error in zeta at each iteration
    niter(1,iter)=iter;
    if iter>=10
        if errorpsi<=ni*nj*0.00000001 && errorzeta<=ni*nj*0.00000001  % Residual Breaking condition
            fprintf("Convergence of Vorticity transport occured at %f iteration\n",iter)
            break
        end
    end
end

%% Solving Pressure Poisson Equation
% Initialisation left bottom corner pressure as 1
p(ni,1)=1;

for iter=1:12000
    
    errorp=0;
    %Left wall
    for row=(ni-1):-1:1
        RHS=-3*zeta(row,1)+4*zeta(row,2)-zeta(row,3);
        old=p(row,1);
        p(row,1)=p(row+1,1)-(RHS/(2*Re));
        errorp=errorp+abs(p(row,1)-old);
    end
    %Bottom wall
    for col=2:nj
        RHS=-3*zeta(ni,col)+4*zeta(ni-1,col)-zeta(ni-2,col);
        old=p(ni,col);
        p(ni,col)=p(ni,col-1)-(RHS/(2*Re));
        errorp=errorp+abs(p(ni,col)-old);
    end
    %Right wall
    for row=(ni-1):-1:1
        RHS=-3*zeta(row,nj)+4*zeta(row,nj-1)-zeta(row,nj-2);
        old=p(row,nj);
        p(row,nj)=p(row+1,nj)-(RHS/(2*Re));
        errorp=errorp+abs(p(row,nj)-old);
    end
    %Pressure Poisson Equation
    for i=2:(nj-1)  % ROW
        for j=2:(ni-1) %COLUMN
            RHS1=(psi(i,j+1)-2*psi(i,j)+psi(i,j-1))*(psi(i-1,j)-2*psi(i,j)+psi(i+1,j));
            RHS2=(psi(i-1,j+1)-psi(i+1,j+1)-psi(i-1,j-1)+psi(i+1,j-1))*(psi(i-1,j+1)-psi(i+1,j+1)-psi(i-1,j-1)+psi(i+1,j-1));
            RHS=RHS1/(power(gridlength,4))-RHS2/(16*power(gridlength,4));
            old=p(i,j);
            p(i,j)=(p(i,j+1)+p(i,j-1)+p(i+1,j)+p(i-1,j)- gridlength*gridlength*(2*RHS))/4;
            errorp=errorp+abs(p(i,j)-old);
        end
    end
    if iter>=10
        if errorp<=ni*nj*0.00000001 %Residual Breaking condition
            fprintf("Convergence of Pressure Poisson occured at %f iteration\n",iter)
            break
        end
    end
end

%% Plots
x = 0:gridlength:length-gridlength;
y = height-gridlength:-gridlength:0;
[X,Y] = meshgrid(x,y);

% Counter of Vorticity values
contour(X,Y,zeta,'ShowText','on')
title('Contour Plot of Vorticity values')
figure;

% Counter of Stream function values
contour(X,Y,psi,'ShowText','on')
title('Contour Plot of Stream function values')
figure;

% Counter of Pressure values
contour(X,Y,p,'ShowText','on')
title('Contour Plot of Pressure values')
figure;

%Plot of Velocity directions
quiver(X,Y,u,v);
title('Plot of Velocity Direction')
xlabel('Length') 
ylabel('Breadth') 
figure;

%Residual Plots
yyaxis left
plot(niter,errorpsimat);
xlabel('No. of Iteration') 
ylabel('Error value of PSI') 
yyaxis right
plot(niter,errorzetamat);
title('Error of psi and zeta vs iteration')
xlabel('No. of Iteration') 
ylabel('Error value of ZETA') 
legend('Psi','Zeta')

            
            
    
   
    
    
            
    
    
            
            
            
    


