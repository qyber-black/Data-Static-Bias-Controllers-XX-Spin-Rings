% script to generate robustness data (currently only logSensitivity) 
% for all controllers in results directory
FILES = dir('results/*.mat')
for k=1:length(FILES)
    robustness = calc_sensitivity(fullfile(FILES(k).folder,FILES(k).name)); 
    save(sprintf('results-robustness/%s-robustness.mat',FILES(k).name(1:end-4)),'robustness'); 
end