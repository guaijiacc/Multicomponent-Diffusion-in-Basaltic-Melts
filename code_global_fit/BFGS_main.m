% this program uses BFGS method to obtain eigenvectors and eigenvalues
% by simultaneously fitting diffusion profiles.

%% Input:
% "****__RealData_8Comp.xlsx":      diffusion data of exps at a temperature ****
% "****__RealError_8Comp.xlsx":     1-sigma errors on oxide concentrations of 
%                                   exps at temperature ****
% "****__RealBoundary_8Comp.xlsx":  oxide concentrations at the left and
%                                   right far field in exps at temperature ****
% "Initial_beta.xlsx":              a l0×7 matrix, with the last 7 rows being 
%                                   eigenvector matrix and the first 3 rows being
%                                   natural log of eigenvalues at three temperatures.

%% Important Outputs:
% S_DF:                             the bulk reduced chi-square;
% trans_beta:                       fitted log of eigenvalues and eigenvector matrix;
% error_trans_beta:                 the 1-sigma errors on "trans_beta";
% data_d(i).RCS:                      the reduced chi-square for the i-th exp;

close all;
clear all;
start_time = clock();

%% ----------------------------define a structure----------------------------
data.temp= 1260; %1260, 1350,1500
data.name = "BS43&44A2";
data.type= 1;% 1:diffusion couple(DC); 2:mineral dissolution(MD)
data.np = 100; % number of points for each exp
data.time = 1000;% exp duration
data.boundary = [];%boundary condition;size:[8,2]
data.x = [];% x-coordinate;[np,1]
data.w = [];% concentration;[np,8]
data.w_calc = []; % predicted concentration;[np,8]; 
data.error=[]; % error of concentration measured by EMPA
data.RCS = 1; % reduced chi-square for each exp

%% ------------------------------initialization---------------------------
T = [1260,1350,1500];%exp temp
% exp duration at different temp
t1 = [1826.0, 1809.6, 1325.5, 1269.2, 950.8, 579.5, 740.0, 577.4, 599.9]; % exp duration at 1260
t2 = [1492.4, 1243.2, 899.0, 974.4, 563.9, 393.6, 522.0, 335.7, 346.3, 441.87]; % exp duration at 1350
t3 = [276.1, 223.6, 247.8, 213.7, 158.5, 157.1, 154.8, 184.9, 221.4]; % exp duration at 1500
t=[t1,t2,t3];

% num of points of each exp
Np1 = [112,121,128,140,97,146,123,128,120]; % number of points of exp at 1260
Np2 = [116,84,100,109,105,101,64,122,123,141];% number of points of exp at 1350
Np3 = [92,124,116,118,120,141,131,141,166]; % number of points of exp at 1500
Np=[Np1,Np2,Np3];

Ndm = [9,10,9]; % number of exp at 1260, 1350 and 1500
threshold = [0,9,19]; % no actual meaning, just for reading data

Nexp = 28; % total number of exp

name1 = ["BS1&2C";"BS3&4C";"BS5&6C";"BS7&8C";"BS9&10C";"BS11&12C";"BS13&14C";"BS17&18C";"BS19&20C"];
name2 = ["BS1&2A";"BS3&4A";"BS5&6A";"BS7&8A";"BS9&10A";"BS11&12A";"BS13&14A";"BS17&18A";"BS19&20A";...
    "BS43&44A2"];
name3 = ["BS1&2B";"BS3&4B";"BS5&6B";"BS7&8B";"BS9&10B";"BS11&12B";"BS13&14B";"BS17&18B";"BS19&20B"];
name = [name1;name2;name3];

%read data
for i=1:3
    T_str = string(T(i));
    realdata = xlsread(strcat(T_str,"_RealData_8Comp.xlsx"));
    realerror = xlsread(strcat(T_str,"_RealError_8Comp.xlsx"));
    realboundary = xlsread(strcat(T_str,"_RealBoundary_8Comp.xlsx"));
    
    for j=1:Ndm(i)
        j1 = j+threshold(i);
        data(j1).temp = T(i);
        data(j1).name = name(j1);
        data(j1).type = 1;
        data(j1).np= Np(j1);
        data(j1).time= t(j1);
        data(j1).boundary = realboundary(1:8,(3*j-2):(3*j-1));
        data(j1).x = realdata(1:Np(j1),12*j-11);
        data(j1).w = realdata(1:Np(j1),(12*j-10):(12*j-3));
        data(j1).error=realerror(:,j);
        
    end
end

% read initial fitting parameters
read_b = xlsread("Initial_beta.xlsx");

%% ---------------------Choose data to fit--------------------------------
% data_d = data;  % Guo data + BS43&44A2
data_d = [data(1:3),data(5:28)];  % Guo data (w/o BS7&8C) + BS43&44A2

[~,m]=size(data_d); % number of exp involved in fitting

Ni = 7; %7 independent components; 

%% ----------------- claculate degree of freedom-------------------------
Npoint = 0;
for i = 1:m
    Npoint = Npoint + data_d(i).np;
end

DF = 8*Npoint - (m+63); % degree of freedom for bulk reduced chi-square;

%% -------------------------initialize error----------------------------
error =[];
for i = 1:m
    Np_i = data_d(i).np;
    error_i = kron(ones(Np_i,1),data_d(i).error);
    error = [error;error_i];
end

Q = read_b(4:10,1:7);% matrix of eigenvector of D
Beta1 = read_b(1,1:7)';% initial log of eigenvalue at 1260 C
Beta2 = read_b(2,1:7)';% initial log of eigenvalue at 1350 C
Beta3 = read_b(3,1:7)';% initial log of eigenvalue at 1500 C


%% ----------------------start of BFGS method-------------------------------
% BFGS method was discovered and published independently by  
% Broyden, Fletcher, Goldfarb, and Shanno in 1970. And it is the most
% famous and widely used quasi-Newton method. Part of variables defined in
% the following part follows LI & FUKUSHIMA, 2001. 

% x, y1, y2 are used for plotting
x = [];
y1 = [];
y2 = [];

n_loop1= 0; % number of outer loops 
n_loop2= 0; % number of inner loops

beta =[ones(m,1);reshape(Q,Ni^2,1);Beta1;Beta2;Beta3];% parameters to be fitted;[m+70,1]

% normalize eigenvectors 
for i =1:7
    beta(m+7*i-6:m+7*i) = beta(m+7*i-6:m+7*i) / norm(beta(m+7*i-6:m+7*i));
end

[n,~] = size(beta);

w0=[]; % initial concentration
for i=1:m
    transform_w = reshape((data_d(i).w)',[],1);
    w0 = [w0;transform_w];
end
w = con(data_d,beta); %calculated concentration
r = (w0-w)./error; % error
S = r'*r; % sum of square error
S_DF = S/DF;   % reduced chi-square
J = calculate_J(data_d,beta); % calculate Jacobian
g = 2*J'*r; %gradient of S (LI & FUKUSHIMA, 2001)
H = eye(m+70); % H matrix (LI & FUKUSHIMA, 2001)

x = [x;n_loop1];
y1 = [y1;S_DF];
y2 = [y2;n_loop2]; 
yyaxis left
plot(x,y1,'-b');
xlabel('num of loop1');
ylabel('reduced chi-square');
yyaxis right
plot(x,y2,'-*r');
ylabel('num of loop2');
drawnow;

max1 = 800;
max2 = 50;
sigma = 0.5; % parameter used in Armijo-type line search (LI & FUKUSHIMA, 2001)
rho = 0.5; % parameter used in Armijo-type line search (LI & FUKUSHIMA, 2001)
n_loop1 = 1;
while n_loop1 < max1 && norm(g)/DF > 1E-6
% when n_loop1 exceeds max1 or ||S_DF|| > 1E-6, iteration ends

    p = -H*g; 
    % the BFGS direction obtained by solving the linear equation: p=-H*g
    step=1; 
    n_loop2 = 0;
    
    % search for proper step using Armijo-type line search.
    while n_loop2 < max2
        beta1 = beta + step*p;
        w1 = con(data_d,beta1);
        r1 = (w0-w1)./error;
        S1 = r1'*r1;
        l= S + sigma*step*g'*p; % parameter used in Armijo-type line search (LI & FUKUSHIMA, 2001)
        
        % if S1<l, the step is OK, otherwise, reduce the step by a factor of rho
        if S1 < l
            break;
        end
        step = rho*step; 
        n_loop2 = n_loop2+1;
    end
    
    if n_loop2==max2
        fprintf("n_loop2 = %d \n",max2);
        break;
    end
    
    s = beta1 - beta; % increase in beta (LI & FUKUSHIMA, 2001)
    % update beta
    beta = beta1;
    % normalize eigenvectors 
    for i =1:7
        beta(m+7*i-6:m+7*i) = beta(m+7*i-6:m+7*i) / norm(beta(m+7*i-6:m+7*i));
    end
    
    % update w,r,S,S_DF,J and y
    w = w1;
    r = r1;
    S = S1;
    S_DF = S/DF;
    J = calculate_J(data_d,beta);
    y = 2*J'*r - g; 
    
    % update H matrix to guarantee the positive definiteness of H
    if norm(g)>1
        if (y'*s)/(s'*s)>= 10^-6*norm(g)^0.01
            H = H + (y'*s + y'*H*y)*(s*s')/(y'*s)^2 - (H*y*s'+s*y'*H)/(y'*s);
        end
    else
        if (y'*s)/(s'*s)>= 10^-6*norm(g)^3
            H = H + (y'*s + y'*H*y)*(s*s')/(y'*s)^2 - (H*y*s'+s*y'*H)/(y'*s);
        end
    end
    
    % update g
    g = 2*J'*r;
 
    x = [x;n_loop1];
    y1 = [y1;S_DF];
    y2 = [y2;n_loop2]; 
    yyaxis left
    plot(x,y1,'-b');
    xlabel('num of loop1');
    ylabel('reduced chi-square');
    yyaxis right
    plot(x,y2,'-*r');
    ylabel('num of loop2');
    drawnow;
    
    n_loop1 = n_loop1 +1;
end

% normalize eigenvectors
Q = reshape(beta(m+1:m+49),7,7);
[~, maxIndices] = max(abs(Q)); % find the indices of the maximum in each column vector of Q
for i =1:7
    Max = beta(m+7*i-7+maxIndices(i));
    beta(m+7*i-6:m+7*i) = beta(m+7*i-6:m+7*i) / norm(beta(m+7*i-6:m+7*i)) * ( Max/abs(Max));
end

%% ------------------------end of BFGS method---------------------------------

transform_w = reshape(w,8,[])';
j=0;

%% store calculated concentration
for i=1:m
    data_d(i).w_calc = transform_w(j+1:j+data_d(i).np,:);
    j = j + data_d(i).np;
end

%% calculate reduced chi-squre for each experiment
for i=1:m
    DF_i = data_d(i).np * 8; % degree of freedom for the i-th exp    
    error_i = kron(ones(data_d(i).np,1),data_d(i).error);
    w0_i = reshape((data_d(i).w)',[],1);
    w_i = reshape((data_d(i).w_calc)',[],1);
    r_i = (w0_i - w_i)./error_i;
    data_d(i).RCS = r_i'*r_i/DF_i;
end

%% calculate eigenvalues
%trans_beta is a 10 by 7 matrix. 
% First row: log of eigenvalues at 1260 
% Second row: log of eigenvalues at 1350 
% Third row: log of eigenvalues at 1500
% The remaining rows: matrix of eigenvector

trans_beta = zeros(10,7);  % change dimensions of beta
trans_beta(1:3,:) = reshape(beta(m+50:m+70) ,7,3)'; % calculate lamda of each exp
trans_beta(4:10,:) = reshape(beta(m+1:m+49),7,7);

%% sequence eigenvalues using bubble sort
for i = 1:6
    for j = 1:7-i
        if trans_beta(1,j)>trans_beta(1,j+1)
            store = trans_beta(:,j+1);
            trans_beta(:,j+1)=trans_beta(:,j);
            trans_beta(:,j)= store;
        end
    end
end

%% calculate D matrix
D1 = trans_beta(4:10,:)*diag(exp(trans_beta(1,:)))/trans_beta(4:10,:);
D2 = trans_beta(4:10,:)*diag(exp(trans_beta(2,:)))/trans_beta(4:10,:);
D3 = trans_beta(4:10,:)*diag(exp(trans_beta(3,:)))/trans_beta(4:10,:);

%% calculate 1-sigma error on trans_beta
beta(m+1:m+49) = reshape(trans_beta(4:10,:),49,1);
beta(m+50:m+70)= reshape(trans_beta(1:3,:)',21,1);
J7 = calculate_J7(data_d,beta);
A = J7'*J7;
inv_A = (diag(inv(A))).^0.5;

error_beta = zeros(m+70,1);
for i = 0:6
    error_beta(m+7*i+1:m+7*i+6)= inv_A(m+6*i+1:m+6*i+6);
end
error_beta(m+50:m+70) = inv_A(m+43:m+63);

error_trans_beta = zeros(10,7);
error_trans_beta(1:3,:) = reshape( error_beta(m+50:m+70) , 7,3)';
error_trans_beta(4:10,:) = reshape( error_beta(m+1:m+49), 7,7);

X0 = beta(1:m); % fitted interface positions
error_X0 = inv_A(1:m); % 1-sigma error on X0

elapsed_time = etime(clock(),start_time); 