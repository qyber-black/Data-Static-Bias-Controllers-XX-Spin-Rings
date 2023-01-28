%% Calcuate log-sensitivy (calc_log_sens)

% SPDX-FileCopyrightText: Copyright (C) 2020-2022 Frank C Langbein <frank@langbein.org>
% SPDX-FileCopyrightText: Copyright (C) 2020-2022 SM Shermer <lw1660@gmail.com>
% SPDX-FileCopyrightText: Copyright (C) 2022 Sean Patrick O'Neil <seanonei@usc.edu>
% SPDX-License-Identifier: CC-BY-SA-4.0  

% This routine takes as input the Hamiltonian for an ring of N spins, a
% controller optimized for maximum fidelity at a given time t or over a 
% time window dt under coherent dynamics, and the time for maximum transfer 
% or time window and calcuates and
% plots the log-sensitivity of the fidelity error to perturbations in the
% Hamiltonian (couplings) or the controls (bias fields). 

% On Input or Used for Execution Computations
% N - size of spin-ring
% type - (1) 'dt' - windowed time transfer given by Info.args.readout
%        (2) 't' - transfer time at t
% out - target spin for transfer (1 -> out)
% obj - Quantum Spin Network (QSN) object for ring of size N that yields
%       the nominal Hamiltonian 
% b - structure of N^2 basis matrices for Bloch transformation 
% D - diagonal matrix of control bias extracted from "data_type_N-out" file
% Hd - controlled Hamiltonian (Hd = obj + D)
% S - structure of 2*N perturbation matrices 
%     S{1} to S{N} are diagonal perturbation matrices 
%     S{N+1} to S{2*N} are perturbations to the couplings 
% t - time of excitation transfer 
% A - N^2 by N^2 matrix of Hd in the Bloch representation 
% lambda - eigenvalues of A
% V - matrix of eigenvectors of A (A = V*lambda*V')
% r_in - Bloch representation of input spin wavefunction 
% r_out - Bloch reprenstation of ouput spin wavefunction 
% c - row vector such that c*exmp(A*t)*r_in provides fidelity err 
% Sb - 2*N strctrue of Bloch transformed perturbation matrices 
% X - N^2 by N^2 matrix of exponentials and 

% On output and saved to
% ../results/log_sens_results/log_sens_data_t/dt-N-out:
% log_sens - array of log-sensitivity of the error filtered for error < 0.1 
%            and calucated for each
%            controller using analytic formula (one row per controller) and
%            ordered by controllers of decreasing fidelity 
%            columns 1 to N are for pertrubations S{1} to S{N}
%            columns N to 2*N are for perturbations S{N+1} to S{2N}
%            column 2*N + 1 has the norm of log-sensitivity for S{k}
% sens - array of sensitivity (partial derivative of the fidelity error to
%        given perturbations) calcuated using the Bloch formulation and
%        ordered by controllers of decreasing fidelity 
%        columns 1 to N are for pertrubations S{1} to S{N}
%        columns N to 2*N are for perturbations S{N+1} to S{2N} 
% log_sens_norm - norm of log_sens for all perturbations for each
%                 controller
% log_sens_bias - norm of log_sens for bias perturbations for each
%                 controller
% log_sens_ham - norm of log_sens for hamiltonian perturbations for each
%                controller                   
% Results - Results field extracted from "data_t/dt-N-out" file
% err - fidelity error taken from the "data_t/dt-N-out" file and ordered
%       by increasing value (decreasing fidelity)

% Figures saved on output:
% log_sens_figure_t/dt_N-out_S' - log-sens vs. fidelity err for each
%                                 perturbation S{k}
% log_sens_composite_t/dt_N-out' - composite log-log snapshot of norm of
% log-sensitivity to Hamilonian perturbations and bias perturbations vs
% err and filtered for fidelity err less than 0.1

clear all; close all; clc;

if ~exist('../results/log_sens_results','dir')
    mkdir('../results/log_sens_results');
end

if ~exist('../figures/log-sensitivity/figures','dir')
    mkdir('../figures/log-sensitivity/figures');
end

if ~exist('../figures/log-sensitivity/composite','dir')
    mkdir('../figures/log-sensitivity/composite');
end

% Outer loop for t vs. dt transfer 
for type = 1:2
    if type == 2
        id = 't';
        start = 2;
    else
        id = 'dt';
        start = 1;
    end

% Loop for ring of each size 
for N = 3:20
    I = eye(N);
    e = arrayfun(@(n) I(:,n),1:N,'UniformOutput',0);
    obj = qsn.QSN('ring',N); % produce QSN object for ring of size N
    b = bloch_basis(N,I);    % produce Bloch basis matrices 
           
    % build perturbation structure matrices 
       for ell = 1:N
           S{ell} = e{ell}*e{ell}';
       end

       for ell = N+1:2*N-1
           S{ell} = e{ell-N}*e{ell-N+1}' + e{ell-N+1}*e{ell-N}';
       end

       S{2*N} = e{1}*e{N}'+e{N}*e{1}';
    
    for out = start:ceil(N/2)
   
    % check for existence of data for transfer and skip if already present  
    if exist(sprintf('../results/log_sens_results/log_sens_data_%s_%d-%d.mat',id,N,out)) ~= 2
       X = zeros(N^2,N^2);   % initialize X for computation of log-sens 
       if out == 1;
           tag1 = sprintf('../data/data_bias_control/data_localisation_%s-%d',id,N);
           tag2 = sprintf('../results/data_bias_control_robustness/data_localisation_%s-%d-robustness',id,N);
       else
           tag1 = sprintf('../data/data_bias_control/data_bias_control_%s-%d-%d',id,N,out);
           tag2 = sprintf('../results/data_bias_control_robustness/data_bias_control_%s-%d-%d-robustness',id,N,out);
       end
   
       load(tag1); % load data files for controller data, transfer time
       load(tag2); 
       
       err = arrayfun(@(x) Results{x}.err,1:2000)';
       time = arrayfun(@(x) Results{x}.time,1:2000)';
       for ell = 1:2000
           controllers(ell,:) = Results{ell}.bias';
       end
       
       % filter controllers for error <0.1 and order 
       Z = [err controllers time];
       Z = sortrows(Z); 
       idx = find(Z(:,1) < 0.1);
       Zfinal = Z(idx,:);
       err = Zfinal(:,1);
       Bias = Zfinal(:,2:N+1);
       time = Zfinal(:,N+2);
  
       % Inner loop for computation of data per transfer 
       M = max(size(err));

       for k = 1:M
       status = sprintf('N=%d target=%d type=%s controller = %d',N,out,id,k);
       disp(status)    

           D = diag(Bias(k,:)); % produce control matrix
           Hd = obj.H+D;              % produce controlled Hamiltonian 
           [A,Sb,lambda,V,r_in,r_out,c] = bloch(N,Hd,S,e,out,b); % Bloch transformation 
           t = time(k,1);       % extract transfer time
           if type == 1
             dt = Info.args.readout(1);
           end

           
           % Inner loop for computation per perturbation  
           for run = 1:2*N            
              Sbar=V'*Sb{run}*V;           % coherent dynamics inv(V) = V'
              
              % calcuate X(t) matrix from analytic formula 
              
              if type == 2 
              
              % Instantaneous transfer calculation

              for m = 1:N^2
                 for n = 1:N^2                    
                    if m == n
                    X(m,n) = Sbar(m,n)*t*exp(lambda(m)*t);
                    else if lambda(m) == lambda(n)
                    X(m,n) = Sbar(m,n)*t*exp(lambda(m)*t);
                    else
                    X(m,n) = (Sbar(m,n)/(lambda(m)-lambda(n)))*(exp(lambda(m)*t)-exp(lambda(n)*t));
                    end
                 end
              end
           end
           
           % set /xi for log-sensitivity computation 
           
           if run <= N
               xi = D(run,run);
           else
               xi = 1;
           end
           
           sens(k,run) = c*V*X*V'*r_in;  % compute sensitivity 
           log_sens(k,run) = abs(xi*c*V*X*V'*r_in)/(err(k,1)); % compute log-sensitivity 
           
     else 
           
          % Time-averaged calculation 
          for m = 1:N^2
           for n = 1:N^2 
            [m n];
           if (lambda(m) == lambda(n)) && (lambda(m) == 0) 
            X(m,n) = Sbar(m,n)*(1/2)*( (t+dt/2)^2 - (t-dt/2)^2 );
           else if (lambda(m) == lambda(n)) && (lambda(m) ~= 0)
            X(m,n) = Sbar(m,n)*(((1/lambda(m))*(t+dt/2)-(1/lambda(m))^2)*exp(lambda(m)*(t+dt/2)) - ((1/lambda(m))*(t-dt/2)-(1/lambda(m))^2)*exp(lambda(m)*(t-dt/2)));
           else if (lambda(m) ~= lambda(n)) && (lambda(m) == 0)
            X(m,n) = (Sbar(m,n)/(lambda(m)-lambda(n)))*( (t+dt/2) - (t-dt/2) - (1/lambda(n))*( exp(lambda(n)*(t + dt/2)) - exp(lambda(n)*(t - dt/2)) ) );
           else if (lambda(m) ~= lambda(n)) && (lambda(n) == 0)
            X(m,n) = (Sbar(m,n)/(lambda(m)-lambda(n)))*( (1/lambda(m))*( exp(lambda(m)*(t + dt/2)) - exp(lambda(m)*(t - dt/2))) - ( (t+dt/2) - (t-dt/2) )  ) ;
           else
            X(m,n) = (Sbar(m,n)/(lambda(m)-lambda(n)))*( (1/lambda(m))*( exp(lambda(m)*(t + dt/2)) - exp(lambda(m)*(t - dt/2))) -  (1/lambda(n))*( exp(lambda(n)*(t + dt/2)) - exp(lambda(n)*(t - dt/2)) ) );
           end
          end
          end
          end
         end
        end

        % set /xi for log-sensitivity computation 
           
        if run <= N
           xi = D(run,run);
        else
           xi = 1;
        end

        sens(k,run) = real(1/dt)*c*V*X*V'*r_in;  % compute sensitivity 
        log_sens(k,run) = abs((1/dt)*(xi*c*V*X*V'*r_in))/(err(k,1)); % compute log-sensitivity 
        end
        end

           log_sens(k,run+1) = norm(log_sens(k,1:run)); % compute norm of log-sensitivity 
       end


log_sens_ham = zeros(M,1);
log_sens_bias = zeros(M,1);
log_sens_norm = zeros(M,1);

% produce data for composite plots 
for k = 1:M
log_sens_ham(k,1) = norm(log_sens(k,N+1:2*N));
log_sens_bias(k,1) = norm(log_sens(k,1:N));
log_sens_norm(k,1) = norm(log_sens(k,2*N+1));
end

% save results 
savetag = sprintf('../results/log_sens_results/log_sens_data_%s_%d-%d',id,N,out);
save(savetag,'log_sens','sens','Results','err','log_sens_ham','log_sens_norm','log_sens_bias');


index = 1:max(size(err));


% plot/save composite figure 
figure;
loglog(err,log_sens_ham,'+',err,log_sens_bias,'*');
grid on;
xlabel('Log(e(t)) [err]');
ylabel('Log(s(\xi_0,t))');
legend('Hamiltonian Perturbations','Bias Perturbations'); 
titletag = sprintf('Log-Sensitivity vs. Fidelity err for %s-%d-%d',id,N,out);
title(titletag);

figtag = sprintf('../figures/log-sensitivity/composite/log_sens_composite_%s_%d-%d',id,N,out);
savefig(figtag);
saveas(gcf,figtag,'png');

% plot/save individual perturbation figures                                         
for ell = 1:2*N

close all;
yyaxis left;
plot(index,log(log_sens(:,ell)));
xlabel('Controller');
ylabel('s(t_f,\xi_0)'); 

yyaxis right;
plot(index,log(err));
ylabel('e(t_f)');
grid on;
set(gcf,'PaperSize',[6 4],'PaperPosition',[0 0 6 4]);
set(findall(gcf,'-property','FontSize'),'FontSize',14);
titletag=sprintf('%s %d-%d Perturbation S=%d',id,N,out,ell);
title(titletag)

figtag = sprintf('../figures/log-sensitivity/figures/log_sens_figure_%s_%d-%d_S=%d',id,N,out,ell);
savefig(figtag);
close all;




end

clear Z;
clear controllers;
clear time;
clear err;
clear log_sens;
clear sens;

end
end
end
end




















