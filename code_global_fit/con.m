% This subroutine calculate concentration profiles for diffusion couple 
% experiments and mineral dissolution experiments with Q (eigenvector matrix) 
% ang log of eigenvalues.

function ww = con(data_d,beta)
    [~,m]=size(data_d);
    Q = reshape(beta(m+1:m+49),7,7);
    ww = [];
    
    for i=1:m
        w1 = data_d(i).boundary(2:8,1); 
        w2 = data_d(i).boundary(2:8,2);
        %for DC, w1 = w_left, w2 = w_right
        %for MD, w1 = w_interface, w2 = w_initial
        
        x0 = beta(i);
        t = data_d(i).time;
        
        if data_d(i).temp == 1260
            Beta = beta(m+50:m+56);
        elseif data_d(i).temp == 1350
            Beta = beta(m+57:m+63);
        else
            Beta = beta(m+64:m+70);
        end
        
        if data_d(i).type==1
            for j=1:data_d(i).np
                x = data_d(i).x(j);
                Y = (x-x0)*exp(-Beta/2)/(sqrt(4*t));
                E = diag( erf(Y) );
                ww1 = (w2 + w1)/2 +Q * E *Q^-1*(w2-w1)/2; %concentration of independent components
                ww2 = 100-sum(ww1,"all"); %concentration of the dependent component,i.e.,SiO2
                ww = [ww;ww2;ww1];
            end
            
        elseif data_d(i).type==2
            for j=1:data_d(i).np
                x = data_d(i).x(j);
                Y = (x-x0)*exp(-Beta/2)/(sqrt(4*t));
                Y0 = -x0*exp(-Beta/2)/(sqrt(4*t));
                E =  diag( erfc(Y)./erfc(Y0) );
                ww1 = w2 +Q * E *Q^-1*(w1-w2); %%concentration of independent components
                ww2 = 100-sum(ww1,"all"); %concentration of the dependent component,i.e.,SiO2
                ww = [ww;ww2;ww1];
            end
        end
    end
end