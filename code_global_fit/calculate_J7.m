% This subroutine calculates the derivates of r(residue) with respect to 
% interface positions, Q(eigenvector matrix) and log of eigenvalues.
 
% J7 differs from J in that J7 doesn't return the derivatives
% of r with respect to the 7th terms of each eigenvector. Therefore, J7 is
% only used to calculate the error while assuming that the errors of the 7th
% terms of each eigenvector are zero.


function J=calculate_J7(data_d,beta)
    [~,m]=size(data_d); %number of experiments included
    Q = reshape(beta(m+1:m+49),7,7); % matrix of eigenvector
    num = 0; 
    n_start = 0;
    n_end = 0;
    JT=[];
    
    for i=1:m
        w1 = data_d(i).boundary(2:8,1);
        w2 = data_d(i).boundary(2:8,2);
        %for DC, w1 = w_left, w2 = w_right
        %for MD, w1 = w_interface, w2 = w_initial
        
        x0 = beta(i);
        t = data_d(i).time;
        error1 = data_d(i).error(2:8);
        error2 = data_d(i).error(1);
        
        if data_d(i).temp == 1260
            num = 0;
            Beta = beta(m+50:m+56);
        elseif data_d(i).temp == 1350
            num = 1;
            Beta = beta(m+57:m+63);
        else
            num = 2;
            Beta = beta(m+64:m+70);
        end
        
        n_start = num*7+1;
        n_end = (num+1)*7;
        
        if data_d(i).type==1   %for diffusion couple
            
            for j=1:data_d(i).np
                x= data_d(i).x(j);
                Y= (x-x0)*exp(-Beta/2)/(sqrt(4*t)); 
                E = diag( erf(Y) ); 
                dE = diag( exp(-Y.*Y).*exp(-Beta/2) ); 
                
                wx0 = -1/sqrt(pi*t)*(Q*dE*Q^-1*(w2-w1)/2)'; % Eq.(18) in "derivation3.pdf"
                % wp is the partial derivative of concentration of 
                % independent components with respect to x0

                wp = kron( (E *Q^-1*(w2-w1)/2), eye(7))-...
                    kron((Q^-1*(w2-w1)/2),(Q*E*Q^-1)'); % Eq.(19) in "derivation3.pdf"
                % wp is the partial derivative of concentration of 
                % independent components with respect to Q
                
                wBeta =-(x-x0)/sqrt(4*pi*t)*diag(reshape(dE,7^2,1))...
                    *kron((Q^-1*(w2-w1)/2),Q'); % Eq.(20) in "derivation3.pdf"
                % wBeta is the partial derivative of concentration of
                % independent components with respect to beta
            
                rx01 =  -wx0/diag(error1);% Eq.(24) in "derivation3.pdf"
                % partial derivative of residue error of independent
                % components with respect with to x0

                rp1 = -wp/diag(error1); % Eq.(25) in "derivation3.pdf"
                % partial derivative of residue error of independent
                % components with respect with to Q
                
                rBeta1= -wBeta/diag(error1); % Eq.(26) in "derivation3.pdf"
                % partial derivative of residue error of independent
                % components with respect with to beta
                
                rx02 = sum(wx0,2)/error2; % Eq.(27) in "derivation3.pdf"
                % partial derivative of residue error of the dependent
                % components with respect with to x0

                rp2 = sum(wp,2)/error2; % Eq.(28) in "derivation3.pdf"
                % partial derivative of residue error of the dependent
                % components with respect with to Q
                
                rBeta2 = sum(wBeta,2)/error2; % Eq.(29) in "derivation3.pdf"
                % partial derivative of residue error of the dependent
                % components with respect with to beta 
                
                rx0 = [rx02,rx01];
                % partial derivative of residue error of all
                % components with respect with to x0

                rp = [rp2,rp1];
                % partial derivative of residue error of all
                % components with respect with to Q
            
                rBeta3=[rBeta2,rBeta1];% 7^2 * (7+1)
                % partial derivative of residue error of all
                % components with respect with to beta

                rX0 = zeros(m,8);
                rX0(i,:) = rx0;
                
                rBeta=zeros(21,8);
                for k=1:8
                    rBeta(n_start:n_end,k)= diag(reshape(rBeta3(:,k),7,7));
                end
                JT =[JT,[rX0;rp;rBeta]];
            end
            
        elseif data_d(i).type==2  % for mineral dissolution
            
            for j=1:data_d(i).np
                x= data_d(i).x(j);
                Y= (x-x0)*exp(-Beta/2)/(sqrt(4*t)); % Eq.(33) in "derivation3.pdf"
                Y0 = -x0*exp(-Beta/2)/(sqrt(4*t)); % Eq.(34) in "derivation3.pdf"
                
                E = erfc(Y); % Eq.(35) in "derivation3.pdf"
                E0 = erfc(Y0); % Eq.(36) in "derivation3.pdf"
                E_ratio = diag (E./E0); % Eq.(37) in "derivation3.pdf"
                 
                dE = exp(-Y.*Y).*exp(-Beta/2); % Eq.(41) in "derivation3.pdf"
                dE0 = exp(-Y0.*Y0).*exp(-Beta/2); % Eq.(43) in "derivation3.pdf"
                
                dE_ratio1 = diag( ( (x-x0)/(sqrt(4*pi*t))*dE.*E0 + x0/(sqrt(4*pi*t))*dE0.*E )./(E0.^2) ); % Eq.(45) in "derivation3.pdf"
                dE_ratio2 = diag( ( 1/(sqrt(pi*t))*dE.*E0 - 1/(sqrt(pi*t))*dE0.*E )./(E0.^2) ); % Eq.(46) in "derivation3.pdf"

                wx0 = (Q*dE_ratio2*Q^-1*(w1-w2))';% Eq.(53) in "derivation3.pdf"
                % wx0 is the partial derivative of concentration of 
                % independent components with respect to x0

                wp = kron( (E_ratio*Q^-1*(w1-w2)), eye(7))-...
                    kron((Q^-1*(w1-w2)),(Q*E_ratio*Q^-1)'); % Eq.(54) in "derivation3.pdf"
                % wp is the partial derivative of concentration of 
                % independent components with respect to Q
                
                wBeta = diag(reshape(dE_ratio1,49,1))...
                    *kron((Q^-1*(w1-w2)),Q'); % Eq.(55) in "derivation3.pdf"
                % wBeta is the partial derivative of concentration of
                % independent components with respect to Beta
                
                rx01 = -wL/diag(error1); % Eq.(59) in "derivation3.pdf"
                % partial derivative of residue error of independent
                % components with respect with to x0

                rp1 = -wp/diag(error1); % Eq.(60) in "derivation3.pdf"
                % partial derivative of residue error of independent
                % components with respect with to Q
                
                rBeta1= -wBeta/diag(error1); % Eq.(61) in "derivation3.pdf"
                % partial derivative of residue error of independent
                % components with respect with to Beta
                
                rx02 = sum(wL,2)/error2; % Eq.(62) in "derivation3.pdf"
                % partial derivative of residue error of the dependent
                % components with respect with to x0
                
                rp2 = sum(wp,2)/error2; % Eq.(63) in "derivation3.pdf"
                % partial derivative of residue error of the dependent
                % components with respect with to Q
                
                rBeta2 = sum(wBeta,2)/error2; % Eq.(64) in "derivation3.pdf"
                % partial derivative of residue error of the dependent
                % components with respect with to Beta 
                
                rx0 = [rx02,rx01];
                % partial derivative of residue error of all
                % components with respect with to x0

                rp = [rp2,rp1];
                % partial derivative of residue error of all
                % components with respect with to Q
                
                rX0 = zeros(m,7+1);
                rX0(i,:) = rx0;

                rBeta3=[rBeta2,rBeta1];%49*8
                % partial derivative of residue error of all
                % components with respect with to Beta 
                
                rBeta=zeros(21,8);
                for k=1:8
                    rBeta(n_start:n_end,k)= diag(reshape(rBeta3(:,k),7,7));
                end
                JT =[JT,[rX0;rp;rBeta]];
            end
        end
    end
    J1 = JT';
    
    [n,~]=size(J1);
    J=zeros(n,m+63);
    for i=0:6
        J(:,m+6*i+1:m+6*i+6)=J1(:,m+7*i+1:m+7*i+6);
    end
    J(:,1:m) = J1(:,1:m);
    J(:,m+43:m+63)=J1(:,m+50:m+70);
end