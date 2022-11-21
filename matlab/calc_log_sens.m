%% Calcuate log-sensitivy (calc_log_sens)

% This routine takes as input the Hamiltonian for an ring of N spins, a
% controller optimized for maximum fidelity at a given time t under 
% coherent dynamics, and the time for maximum transfer and calcuates and
% plots the log-sensitivity of the fidelity error to perturbations in the
% Hamiltonian (couplings) or the controls (bias fields). 

% On Input or Used for Execution Computations
% N - size of spin-ring
% type - (1) 'dt' - time transfer of window of 0.1 about t
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
% c - row vector such that c*exmp(A*t)*r_in provides fidelity error 
% Sb - 2*N strctrue of Bloch transformed perturbation matrices 
% X - N^2 by N^2 matrix of exponentials and 

% On output and saved to
% log-sensitivity-new/log_sens_results/log_sens_data_t/dt-N-out:
% log_sens - array of log-sensitivity of the error calucated for each
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
% log_sens_filtered - copy of log-sens array but with all data for
%                     controllers producing a fielity error greater than 0.1 removed 
% log_sens_short - structre with three arrays that summarize the filtered
%                  log-sensitiviy data 
%                  .ham - norm of log_sens_filtered for Hamiltonian
%                  perturbations for each controller 
%                  .bias - norm of log_sens_filtered for bias perturbations
%                  .norm - norm of log_sens_filtered for all perturbations 
% results - results field extracted from "data_t/dt-N-out" file
% error1 - fidelity error taken from the "data_t/dt-N-out" file and ordered
%          by increasing value (decreasing fidelity)
% error2 - fidelity error calculated by computing error at transfer time t
%          in the Bloch formulation 
% old_sens - sensitivity extracted from the "data_t/dt-N-out" file and
%            ordered by decreasing fidelity 
 

% Figures saved on output:
% log_sens_figure_t/dt_N-out_S' - log-sens vs. fidelity error for each
%                                 perturbation S{k}
% log_sens_composite_t/dt_N-out' - composite log-log snapshot of norm of
% log-sensitivity to Hamilonian perturbations and bias perturbations vs
% error and filtered for fidelity error less than 0.1



clear all; close all; clc;

if ~exist('../log-sensitivity-new/log_sens_results','dir')
    mkdir('../log-sensitivity-new/log_sens_results');
end

if ~exist('../log-sensitivity-new/figures','dir')
    mkdir('../log-sensitivity-new/figures');
end

if ~exist('../log-sensitivity-new/figures/composite','dir')
    mkdir('../log-sensitivity-new/figures/composite');
end

% Outer loop for t vs. dt transfer 
for type = 2:2
    if type == 2
        id = 't';
        start = 2;
    else
        id = 'dt';
        start = 1;
    end

% Loop for ring of each size 
for N = 11:20
    I = eye(N);
    e = arrayfun(@(n) I(:,n),1:N,'UniformOutput',0);
    obj = qsn.QSN('ring',N); % produce QSN object for ring of size N
    b = bloch_basis(N,I);    % produce Bloch basis matrices 
    for out = start:ceil(N/2)
       X = zeros(N^2,N^2);   % initialize X for computation of log-sens 
       tag1 = sprintf('data/data_%s-%d-%d',id,N,out);
       if out == 1;
           tag2 = sprintf('results-robustness/data_localisation_%s-%d-robustness',id,N);
       else
           tag2 = sprintf('results-robustness/data_bias_control_%s-%d-%d-robustness',id,N,out);
       end
       load(tag1); % load data files to extract controllers and transfer time
       load(tag2); % load old robustness files for data comparison 
       
       % build perturbation structure matrices 
       for ell = 1:N
           S{ell} = e{ell}*e{ell}';
       end

       for ell = N+1:2*N-1
           S{ell} = e{ell-N}*e{ell-N+1}' + e{ell-N+1}*e{ell-N}';
       end

       S{2*N} = e{1}*e{N}'+e{N}*e{1}';
        
       % Inner loop for computation of data per transfer 
       M = max(size(results));
       status = sprintf('N=%d target=%d type=%s',N,out,id);
       disp(status);
       for k = 1:M
           D = diag(results{k}.bias); % produce control matrix
           Hd = obj.H+D;              % produce controlled Hamiltonian 
           [A,Sb,lambda,V,r_in,r_out,c] = bloch(N,Hd,S,e,out,b); % Bloch transformation 
           t = results{k}.time;       % extract transfer time
           error1(k,1) = results{k}.err;    % extract original fidelity error 
           error2(k,1) = c*expm(A*t)*r_in;  % compute updated fidelity error 
           
           % Inner loop for computation per perturbation 
           for run = 1:2*N            
              Sbar=V'*Sb{run}*V;           % coherent dynamics inv(V) = V'
              % calcuate X(t) matrix from analytic formula 
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
           log_sens(k,run) = abs(xi*c*V*X*V'*r_in)/(error1(k,1)); % compute log-sensitivity 
           end
           log_sens(k,run+1) = norm(log_sens(k,1:run)); % compute norm of log-sensitivity 
end

% order log-sensitivity and sensitivity by decreasing fidelity 
Z = [error1 error2 log_sens sens];
Z = sortrows(Z); 

error1 = Z(:,1);
error2 = Z(:,2);
log_sens = Z(:,3:2*N+3);
sens = Z(:,2*N+4:4*N+3);

% order old sensitivity by decreasing fidelity 
Z3 = [robustness.error' cell2mat(robustness.sensitivity)'];
Z3 = sortrows(Z3);
old_sens = Z3(:,2:2*N+1); 


% filter log-sensitivity for high fidelity controllers 
idx = find(error1 < 0.1);
Z2 = Z(idx,:);
log_sens_filtered = Z2(:,3:2*N+2);
error_filtered = Z2(:,1);
M2 = max(size(error_filtered));

log_sens_short.ham = zeros(M2,1);
log_sens_short.bias = zeros(M2,1);
log_sens_short.norm = zeros(M2,1);

% produce data for composite plots 
for k = 1:M2
log_sens_short.ham(k,1) = norm(log_sens_filtered(k,N+1:2*N));
log_sens_short.bias(k,1) = norm(log_sens_filtered(k,1:N));
log_sens_short.norm(k,1) = norm(log_sens_filtered(k,1:2*N));
end

% save results 
savetag = sprintf('../log-sensitivity-new/log_sens_results/log_sens_data_%s_%d-%d',id,N,out);
save(savetag,'log_sens','sens','results','error1','error2','log_sens_filtered','log_sens_short','old_sens');

index = 1:max(size(error1));

% plot/save composite figure 
figure;
loglog(error_filtered,log_sens_short.ham,'+',error_filtered,log_sens_short.bias,'*');
grid on;
xlabel('Log(e(t_f)) [error at readout time t_f]');
ylabel('Log(s(\xi_0,t_f))');
legend('Hamiltonian Perturbations','Bias Perturbations'); 
titletag = sprintf('Log-Sensitivity vs. Fidelity Error for %s-%d-%d',id,N,out);
title(titletag);

figtag = sprintf('../log-sensitivity-new/figures/composite/log_sens_composite_%s_%d-%d',id,N,out);
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
plot(index,log(error1));
ylabel('e(t_f)');
grid on;
set(gcf,'PaperSize',[6 4],'PaperPosition',[0 0 6 4]);
set(findall(gcf,'-property','FontSize'),'FontSize',14);
titletag=sprintf('%s %d-%d Perturbation S=%d',id,N,out,ell);
title(titletag)

figtag = sprintf('../log-sensitivity-new/figures/log_sens_figure_%s_%d-%d_S=%d',id,N,out,ell);
savefig(figtag);




end
end
end
end




















