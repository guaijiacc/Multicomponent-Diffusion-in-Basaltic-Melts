%solving a diffusion eq for a diffusion couple in which D depends on conc
% using the explicit method; D = D1*exp(a*C);
% dC/dt=d(D1*exp(a*C)*dC/dx)/dx; D1 is D at C=0
% x = 0 is the Matano interface
% initial condition C=Cmin at x<0; C=Cmax at x>0
% ratio=D1*(dt)/(dx)^2
% fit parameters: D1; a; x0 (effective binary diffusion)
% for different profiles, it is necessary to change the following:
%  (1) Input and output data file name, and ylabel; 
%  (2) Cmin, Cmax, dt, ttotal (numt*dt=duration);
%  (3) xmin, xmax, dx (make dx = 10 µm or less);
%  (4) initial guesses of D1 (D1 is minimum D), a, x0
% Inputs need to be checked if the line ends with "% text"

tic
close all; 
clear all;

xmin=-1500; xmax=1500; dx=5; %range of x in micrometers
numx=1+(xmax-xmin)/dx; 
numx1=1+(numx-1)*(0-xmin)/(xmax-xmin); %space step at x=0 (interface)
x = xmin:dx:xmax;  
dt=.01; ttotal=740; % time step and total time (exp duration) in s
numt=ttotal/dt;  
t = 0:dt:ttotal;
Np=5; % number of fitting parameters

%% For Z6
K = 6;
[Zi,text,all] = xlsread('BS13&14C_Z6.xlsx'); % data (x,C) to be fit
Par = [2.749, 0.183, 16.685, 9.538, 14.945]'; % 5 fitting parameters for Z6

D1 = Par(1);
a = Par(2);
x0 = Par(3);
Cmin = Par(4);
Cmax = Par(5);

[m1,n1]=size(Zi);
Par_err=zeros(Np,1); V1=zeros(2,1);
Fp=zeros(m1,1); Fd=zeros(m1,Np);Temp1=0; SStot=0;
for i=1:m1
    Temp1=Temp1+Zi(i,2);
end
Temp1=Temp1/m1;
for i=1:m1
    SStot=SStot+(Zi(i,2)-Temp1)^2;
end

format long;

delta = 1E-2; % adjust the step as needed!

for k=1:20
    [X1] = Solver3(D1,a,Cmin,Cmax,x,dx,dt,numx,numx1,numt);
    % Interpolate to obtain K2O concentration at each analyzed point
    xdata=Zi(:,1)-x0;
    ycalc=interp1(X1(:,1),X1(:,2),xdata);
    Fp=Zi(:,2)-ycalc; Chisq=Fp'*Fp; r_squared=1-Chisq/SStot;
    % Fp1=Fp./Zi(:,3); Chisq1=Fp1'*Fp1;
    disp(Par'); V1(1)=sqrt(Chisq/(m1-Np)); V1(2)=r_squared;
    disp(V1'); 

    figure(1);
    hold on;
    xdata1=X1(:,1)+x0;
    plot(xdata1,X1(:,2));
    drawnow;
    title('Diffusion couple profile; D=D1*exp(a*C)');
    xlabel('x(micrometer)');
    ylabel('Z6');
    plot(Zi(:,1),Zi(:,2),'o');
    drawnow;

    Par_a=Par; Par_b=Par;
    %   Par_a(1)=Par(1,1)*0.99; Par_b(1)=Par(1,1)*1.01; dP=Par_b(1)-Par_a(1);
    %   xcalc_a=X1(:,1)*sqrt(0.99); xcalc_b=X1(:,1)*sqrt(1.01);
    Par_a(1)=Par(1,1)*(1-delta); Par_b(1)=Par(1,1)*(1+delta); dP=Par_b(1)-Par_a(1);
    xcalc_a=X1(:,1)*sqrt((1-delta)); xcalc_b=X1(:,1)*sqrt((1+delta));
    ycalc_a=interp1(xcalc_a,X1(:,2),xdata);
    ycalc_b=interp1(xcalc_b,X1(:,2),xdata);
    Fp_a=Zi(:,2) - ycalc_a; % Fp1_a=Fp_a./Zi(:,3); 
    Fp_b=Zi(:,2) - ycalc_b; % Fp1_b=Fp_b./Zi(:,3); 
    Fd(:,1)=(Fp_b-Fp_a)/dP;

    Par_a=Par; Par_b=Par;
    %   Par_a(2)=Par(2,1)*0.99; Par_b(2)=Par(2,1)*1.01; dP=Par_b(2)-Par_a(2);
    Par_a(2)=Par(2,1)*(1-delta); Par_b(2)=Par(2,1)*(1+delta); dP=Par_b(2)-Par_a(2);
    [X1a] = Solver3(Par_a(1),Par_a(2),Par_a(4),Par_a(5),x,dx,dt,numx,numx1,numt);
    ycalc_a=interp1(X1a(:,1),X1a(:,2),xdata);
    [X1b] = Solver3(Par_b(1),Par_b(2),Par_b(4),Par_b(5),x,dx,dt,numx,numx1,numt);
    ycalc_b=interp1(X1b(:,1),X1b(:,2),xdata);
    Fp_a=Zi(:,2) - ycalc_a; % Fp1_a=Fp_a./Zi(:,3); 
    Fp_b=Zi(:,2) - ycalc_b; % Fp1_b=Fp_b./Zi(:,3); 
    Fd(:,2)=(Fp_b-Fp_a)/dP;

    Par_a=Par; Par_b=Par;
    %   Par_a(3)=Par(3,1)-0.02; Par_b(3)=Par(3,1)+0.02; dP=Par_b(3)-Par_a(3);
    Par_a(3)=Par(3,1)*(1-delta); Par_b(3)=Par(3,1)+(1+delta); dP=Par_b(3)-Par_a(3);
    xdata_a=Zi(:,1)-Par_a(3);
    ycalc_a=interp1(X1(:,1),X1(:,2),xdata_a);
    xdata_b=Zi(:,1)-Par_b(3);
    ycalc_b=interp1(X1(:,1),X1(:,2),xdata_b);
    Fp_a=Zi(:,2) - ycalc_a; % Fp1_a=Fp_a./Zi(:,3); 
    Fp_b=Zi(:,2) - ycalc_b; % Fp1_b=Fp_b./Zi(:,3); 
    Fd(:,3)=(Fp_b-Fp_a)/dP;

    Par_a=Par; Par_b=Par;
    %   Par_a(4)=Par(4,1)*0.99; Par_b(4)=Par(4,1)*1.01; dP=Par_b(4)-Par_a(4);
    Par_a(4)=Par(4,1)*(1-delta); Par_b(4)=Par(4,1)*(1+delta); dP=Par_b(4)-Par_a(4);
    [X1a] = Solver3(Par_a(1),Par_a(2),Par_a(4),Par_a(5),x,dx,dt,numx,numx1,numt);
    ycalc_a=interp1(X1a(:,1),X1a(:,2),xdata);
    [X1b] = Solver3(Par_b(1),Par_b(2),Par_b(4),Par_b(5),x,dx,dt,numx,numx1,numt);
    ycalc_b=interp1(X1b(:,1),X1b(:,2),xdata);
    Fp_a=Zi(:,2) - ycalc_a; % Fp1_a=Fp_a./Zi(:,3); 
    Fp_b=Zi(:,2) - ycalc_b; % Fp1_b=Fp_b./Zi(:,3); 
    Fd(:,4)=(Fp_b-Fp_a)/dP;

    Par_a=Par; Par_b=Par;
    %   Par_a(5)=Par(5,1)*0.99; Par_b(5)=Par(5,1)*1.01; dP=Par_b(5)-Par_a(5);
    Par_a(5)=Par(5,1)*(1-delta); Par_b(5)=Par(5,1)*(1+delta); dP=Par_b(5)-Par_a(5);
    [X1a] = Solver3(Par_a(1),Par_a(2),Par_a(4),Par_a(5),x,dx,dt,numx,numx1,numt);
    ycalc_a=interp1(X1a(:,1),X1a(:,2),xdata);
    [X1b] = Solver3(Par_b(1),Par_b(2),Par_b(4),Par_b(5),x,dx,dt,numx,numx1,numt);
    ycalc_b=interp1(X1b(:,1),X1b(:,2),xdata);
    Fp_a=Zi(:,2) - ycalc_a; % Fp1_a=Fp_a./Zi(:,3); 
    Fp_b=Zi(:,2) - ycalc_b; % Fp1_b=Fp_b./Zi(:,3); 
    Fd(:,5)=(Fp_b-Fp_a)/dP;

    % Newtonian method
    dPar=inv(Fd'*Fd)*Fd'*Fp;
    %   Par=Par-0.5*dPar; 
    Par=Par-dPar; 
    D1=Par(1,1); a=Par(2,1);x0=Par(3,1); Cmin=Par(4,1); Cmax=Par(5,1);
end

%Calculation of errors (weight on measured conc is identical; no error on x)
%Clifford 1973 book Multivariate Error Analysis p.79
Sigma=sqrt((Fp'*Fp)/(m1-Np)); P=inv(Fd'*Fd);
for i=1:Np
    Par_err(i)=Sigma*sqrt(P(i,i));
end

format long;
disp('r^2='); disp(r_squared);
disp('mean error on Zi'); disp(sqrt(Chisq/(m1-Np)));
disp('solution (D1, a, x0, Cmin, Cmax)'); disp(Par');
disp('error'); disp(Par_err');

%% write fitted Zi into file
variable_name = ["x(μm)","Z" + string(K) + "_fit"];
filename = "fit_Z" + string(K) + "_dependent.xlsx";
T1 = table(xdata1,X1(:,2));
T = splitvars(T1);
T.Properties.VariableNames = variable_name;
writetable(T,filename,"AutoFitWidth",false);

%% BS13&14C, Z6, fitting result
% r^2: 0.999410234367520
% mean error on Zi: 0.056885874157884
% [D1, a, x0, Cmin, Cmax]:
% [2.751469942215767   0.183224938734134  16.682156769746271   9.538393403144301  14.944686184932726]
% error:[0.492298537873168   0.014693742035484   2.010889119356218   0.010365565269915   0.012565047090866]

toc