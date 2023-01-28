%% Analyze log-sensitivity for -t controllers (analyze_log_sens_t)

% SPDX-FileCopyrightText: Copyright (C) 2020-2022 Frank C Langbein <frank@langbein.org>
% SPDX-FileCopyrightText: Copyright (C) 2020-2022 SM Shermer <lw1660@gmail.com>
% SPDX-FileCopyrightText: Copyright (C) 2022 Sean Patrick O'Neil <seanonei@usc.edu>
% SPDX-License-Identifier: CC-BY-SA-4.0  

% This script takes as input the controller data in the
% data_bias_control_t_N-T files and the robustness data from the
% calc_sensitivity.m routine saved in the
% data_bias_control_dt_N-T-robustness.m files and computes the correlation
% between the error and log-sensitivity along with figures depicting the
% relation. 

% On Input
% N - size of spin-ring
% target - target spin for transfer (1 -> out)
 

% On output
% Corrleation_Data.xlsx - spreadsheet of correlation data organized in the
% following columns:
%   1 - Kendall tau for error and norm of log-sensitivity 
%   2 - Kendall tau test statistic for (1)
%   3 - p-value for test statistic in (3) 
%   4 - Kendall tau for error and norm of log-sensitivity to control uncertainty  
%   5 - Kendall tau test statistic for (4)
%   6 - p-value for test statistic in (5) 
%   7 - Kendall tau for error and norm of log-sensitivity to hamiltonian uncertainty 
%   8 - Kendall tau test statistic for (7)
%   9 - p-value for test statistic in (8) 
%   10 - Pearson r for log(error) and log(norm(log-sensitivity)) 
%   11 - Pearson r test statistic for (10)
%   12 - p-value for test statistic in (11) 
%   13 - Pearson r for log(error) and log(norm(log-sensitivity)) to control uncertainty  
%   14 - Pearson r test statistic for (13)
%   15 - p-value for test statistic in (14) 
%   16 - Pearson r for log(error) and log(norm(log-sensitivity)) to hamiltonian uncertainty 
%   17 - Pearson r test statistic for (16)
%   18 - p-value for test statistic in (17) 
%
% /figures/log-sensitivity/comparison/log_sens_composite_t_%d-%d_comparison
%   figures showing a log-log plot of the norm of the log-sensitivity and norm of 
%   log-sensitivity to hamiltonian and controller uncertainty versus error
%   providing a visual depiction of the correlation data 

clear all; close all; clc;

if ~exist('../figures/log-sensitivity/comparison','dir')
    mkdir('../figures/log-sensitivity/comparison');
end

if ~exist('../results/log_sens_results','dir')
    mkdir('../figures/log_sens_results');
end

count = 0;

for N = 3:20
for target = 2:ceil(N/2)
    count = count+1;
        tag1 = sprintf('../results/data_bias_control_robustness/data_bias_control_t-%d-%d-robustness',N,target);
        file = sprintf('../data/data_bias_control/data_bias_control_t-%d-%d',N,target);
    load(tag1); %load robustness data 
    load(file); %load controller data 
    I = eye(N);
    e = arrayfun( @(x) I(:,x), 1:N, 'UniformOutput',false);
    error = robustness.error';
    
    % extract log-sensitivity 
    for k = 1:length(error)
        log_sens(k,:) = robustness.log_sensitivity{k};
    end

    for k = 1:length(error)
        for l = 1:N
                log_sens(k,l) = log_sens(k,l)*Results{k}.bias(l); %scale bias uncertainty by nominal value of bias field 
        end
        log_sens(k,2*N+1) = norm(log_sens(k,:));
        log_sens_hamil(k,1) = norm(log_sens(k,N+1:2*N));
        log_sens_control(k,1) = norm(log_sens(k,1:N));
    end

    L = [error log_sens log_sens_control log_sens_hamil];
    L2 = sortrows(L,'ascend');
    idx = find(L2(:,1)>0.1);
    L2(idx,:) = [];
    error = L2(:,1);
    log_sens = L2(:,2:2*N+2);
    log_sens_hamil = L2(:,2*N+4);
    log_sens_control = L2(:,2*N+3);
    
    % Calculate Kendall tau statistics
    correl(count,1) = corr(error,log_sens(:,2*N+1),'type','kendall');
    correl(count,4) = corr(error,log_sens_control,'type','kendall');
    correl(count,7) = corr(error,log_sens_hamil,'type','kendall');
    n = length(error);
    sigma_k = sqrt(2*(2*n+5)/(9*n*(n-1)));
    correl(count,2) = correl(count,1)/sigma_k;
    correl(count,5) = correl(count,4)/sigma_k;
    correl(count,8) = correl(count,7)/sigma_k;
    
    if correl(count,2) > 0
        correl(count,3) = (1-normcdf(correl(count,2)));
    else
        correl(count,3) = normcdf(correl(count,2));
    end
    
    if correl(count,5) > 0
        correl(count,6) = (1-normcdf(correl(count,5)));
    else
        correl(count,6) = normcdf(correl(count,5));
    end

    if correl(count,8) > 0
        correl(count,9) = (1-normcdf(correl(count,8)));
    else
        correl(count,9) = normcdf(correl(count,8));
    end

    % Calculate Pearson r statistics 
    correl(count,10) = corr(log(error),log(log_sens(:,2*N+1)),'type','pearson');
    correl(count,13) = corr(log(error),log(log_sens_control),'type','pearson');
    correl(count,16) = corr(log(error),log(log_sens_hamil),'type','pearson');
    n = length(error);
    sigma_p = sqrt(n-2);
    correl(count,11) = correl(count,10)*sigma_p/sqrt(1 - correl(count,10)^2);
    correl(count,14) = correl(count,13)*sigma_p/sqrt(1 - correl(count,13)^2);
    correl(count,17) = correl(count,16)*sigma_p/sqrt(1 - correl(count,16)^2);
    
    if correl(count,11) > 0
        correl(count,12) = (1-tcdf(correl(count,11),n-2));
    else
        correl(count,12) = tcdf(correl(count,11),n-2);
    end
    
    if correl(count,14) > 0
        correl(count,15) = (1-tcdf(correl(count,14),n-2));
    else
        correl(count,15) = tcdf(correl(count,14),n-2);
    end

    if correl(count,17) > 0
        correl(count,18) = (1-tcdf(correl(count,17),n-2));
    else
        correl(count,18) = tcdf(correl(count,17),n-2);
    end
    
        correl(count,19) = n;

    rowtag = sprintf('N=%d out=%d',N,target);
    rowname{count,1} = rowtag;
   
    
    % plot/save composite figure 
    index = length(error);
    figure;
    loglog(error,log_sens(:,2*N+1),'+',error,log_sens_control(:,1),'*',error,log_sens_hamil,'.');
    grid on;
    xlabel('Log(err) [error at readout time]');
    ylabel('Log(s(\xi_0,t_f))');
    legend('Norm of log-sensitivity','Log-sensitivity - Control Perturbations','Log-sensitivity - Hamiltonian Perturbations'); 
    titletag = sprintf('Log-Sensitivity vs. Fidelity Error for t-%d-%d',N,target);
    title(titletag);

    figtag = sprintf('../figures/log-sensitivity/comparison/log_sens_composite_t_%d-%d_comparison',N,target);
    savefig(figtag);
    saveas(gcf,figtag,'png');
    

    clear log_sens;
    close;
end
end

headings = {'tau - norm(log_sens) v. err','Test Stat1','p-value1','tau - norm(log_sens_control) vs. err','Test Stat2','p-value2','tau - norm(log_sens_hamil) vs. err','Test Stat3','p-value3','Pearson - log(norm(log_sens)) v. log(err)','Test-Stat4','p-value4','Pearson - log(norm(log_sens_control)) vs. log(err)','Test Stat5','p-value5','Pearson - log(norm(log_sens_hamil)) vs. log(err)','Test Stat6','p-value6','n'}

corr_data = array2table(correl);
corr_data.Properties.RowNames = rowname;
corr_data.Properties.VariableNames = headings;
writetable(corr_data,'../results/log_sens_results/Correlation_Data.xlsx','sheet','t_Controllers','WriteRowNames',true);
