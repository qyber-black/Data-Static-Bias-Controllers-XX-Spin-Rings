%% Analyze log-sensitivity for dt controllers (analyze_log_sens_dt)

% SPDX-FileCopyrightText: Copyright (C) 2020-2022 Frank C Langbein <frank@langbein.org>
% SPDX-FileCopyrightText: Copyright (C) 2020-2022 SM Shermer <lw1660@gmail.com>
% SPDX-FileCopyrightText: Copyright (C) 2022 Sean Patrick O'Neil <seanonei@usc.edu>
% SPDX-License-Identifier: CC-BY-SA-4.0  

% This script takes as input the controller data in the
% data_bias_control_dt_N-T files and the robustness data from the
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
% /figures/log-sensitivity/comparison/log_sens_composite_dt_%d-%d_comparison
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
for N = 3:20;
for target = 1:ceil(N/2)
    count = count+1;
    if target ~=1
        tag1 =sprintf('../results/data_bias_control_robustness/data_bias_control_dt-%d-%d-robustness',N,target);
        file = sprintf('../data/data_bias_control/data_bias_control_dt-%d-%d',N,target);
    else 
        tag1 = sprintf('../results/data_bias_control_robustness/data_localisation_dt-%d-robustness',N);
        file = sprintf('../data/data_bias_control/data_localisation_dt-%d',N);
    end

    load(tag1);
    load(file);
    I = eye(N);
    e = arrayfun( @(x) I(:,x), 1:N, 'UniformOutput',false);

    
    error_dt = robustness.error';

    for k = 1:length(error_dt)
        log_sens_dt(k,:) = robustness.log_sensitivity{k};
    end

    for k = 1:length(error_dt)
        for l = 1:N
                log_sens_dt(k,l) = log_sens_dt(k,l)*Results{k}.bias(l);
        end
        log_sens_dt(k,:) = abs(log_sens_dt(k,:));
        log_sens_dt(k,2*N+1) = norm(log_sens_dt(k,:));
        log_sens_dt_hamil(k,1) = norm(log_sens_dt(k,N+1:2*N));
        log_sens_dt_control(k,1) = norm(log_sens_dt(k,1:N));

    end

    Z = [error_dt log_sens_dt log_sens_dt_control log_sens_dt_hamil];
    Z2 = sortrows(Z,'ascend');
    idx = find(Z2(:,1) > 0.1)
    Z2(idx,:) = [];
    error_dt = Z2(:,1);
    log_sens_dt = Z2(:,2:2*N+2);
    log_sens_dt_hamil = Z2(:,2*N+4);
    log_sens_dt_control = Z2(:,2*N+3);


        % Calculate Kendall tau statistics
    corr_dt(count,1) = corr(error_dt,log_sens_dt(:,2*N+1),'type','kendall');
    corr_dt(count,4) = corr(error_dt,log_sens_dt_control,'type','kendall');
    corr_dt(count,7) = corr(error_dt,log_sens_dt_hamil,'type','kendall');
    n = length(error_dt);
    sigma_k = sqrt(2*(2*n+5)/(9*n*(n-1)));
    corr_dt(count,2) = corr_dt(count,1)/sigma_k;
    corr_dt(count,5) = corr_dt(count,4)/sigma_k;
    corr_dt(count,8) = corr_dt(count,7)/sigma_k;
    
    if corr_dt(count,2) > 0
        corr_dt(count,3) = 1*(1-normcdf(corr_dt(count,2)));
    else
        corr_dt(count,3) = 1*normcdf(corr_dt(count,2));
    end
    
    if corr_dt(count,5) > 0
        corr_dt(count,6) = 1*(1-normcdf(corr_dt(count,5)));
    else
        corr_dt(count,6) = 1*normcdf(corr_dt(count,5));
    end

    if corr_dt(count,8) > 0
        corr_dt(count,9) = 1*(1-normcdf(corr_dt(count,8)));
    else
        corr_dt(count,9) = 1*normcdf(corr_dt(count,8));
    end

    % Calculate Pearson r statistics 
    corr_dt(count,10) = corr(log10(error_dt),log10(log_sens_dt(:,2*N+1)),'type','pearson');
    corr_dt(count,13) = corr(log10(error_dt),log10(log_sens_dt_control),'type','pearson');
    corr_dt(count,16) = corr(log10(error_dt),log10(log_sens_dt_hamil),'type','pearson');
    n = length(error_dt);
    sigma_p = sqrt(n-2);
    corr_dt(count,11) = corr_dt(count,10)*sigma_p/sqrt(1 - corr_dt(count,10)^2);
    corr_dt(count,14) = corr_dt(count,13)*sigma_p/sqrt(1 - corr_dt(count,13)^2);
    corr_dt(count,17) = corr_dt(count,16)*sigma_p/sqrt(1 - corr_dt(count,16)^2);
    
    if corr_dt(count,11) > 0
        corr_dt(count,12) = 1*(1-tcdf(corr_dt(count,11),n-2));
    else
        corr_dt(count,12) = 1*tcdf(corr_dt(count,11),n-2);
    end
    
    if corr_dt(count,14) > 0
        corr_dt(count,15) = 1*(1-tcdf(corr_dt(count,14),n-2));
    else
        corr_dt(count,15) = 1*tcdf(corr_dt(count,14),n-2);
    end

    if corr_dt(count,17) > 0
        corr_dt(count,18) = 1*(1-tcdf(corr_dt(count,17),n-2));
    else
        corr_dt(count,18) = 1*tcdf(corr_dt(count,17),n-2);
    end

    corr_dt(count,19) = n;
    


    rowtag = sprintf('N=%d out=%d Transfer',N,target);
    rowname{count,1} = rowtag;
    
    % plot/save composite figure 
    index = length(error_dt);
    figure;
    loglog(error_dt,log_sens_dt(:,2*N+1),'+',error_dt,log_sens_dt_control(:,1),'*',error_dt,log_sens_dt_hamil,'.');
    grid on;
    xlabel('Log(err) [error over readout window]');
    ylabel('Log(s(\xi_0,t_f))');
    legend('Norm of log-sensitivity','Log-sensitivity - Control Perturbations','Log-sensitivity - Hamiltonian Perturbations'); 
    titletag = sprintf('Log-Sensitivity vs. Fidelity Error for dt-%d-%d',N,target);
    title(titletag);

    figtag = sprintf('../figures/log-sensitivity/comparison/log_sens_composite_dt_%d-%d_comparison',N,target);
    savefig(figtag);
    saveas(gcf,figtag,'png');
    
    clear log_sens_dt;
    close;
end
end

headings = {'tau - norm(log_sens) v. err','Test Stat1','p-value1','tau - norm(log_sens_control) vs. err','Test Stat2','p-value2','tau - norm(log_sens_hamil) vs. err','Test Stat3','p-value3','Pearson - log(norm(log_sens)) v. log(err)','Test-Stat4','p-value4','Pearson - log(norm(log_sens_control)) vs. log(err)','Test Stat5','p-value5','Pearson - log(norm(log_sens_hamil)) vs. log(err)','Test Stat6','p-value6','n'}

corr_data_dt = array2table(corr_dt);
corr_data_dt.Properties.RowNames = rowname;
corr_data_dt.Properties.VariableNames = headings;
writetable(corr_data_dt,'../results/log_sens_results/Correlation_Data.xlsx','sheet','dt_Controllers','WriteRowNames',true);
