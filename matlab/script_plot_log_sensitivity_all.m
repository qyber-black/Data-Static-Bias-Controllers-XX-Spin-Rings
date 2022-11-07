% script to plot log sensitivity data
FILES = dir('results-robustness/*.mat')
for k=1:length(FILES)
    load(fullfile(FILES(k).folder,FILES(k).name)); 
    figure(1), clf, 
    plot_logSens_vs_error(robustness), drawnow
    title(FILES(k).name)
    savefig(1,sprintf('../figures/log-sensitivity/%s',FILES(k).name(1:end-4))); 
    saveas(1,sprintf('../figures/log-sensitivity/%s',FILES(k).name(1:end-4)),'png'); 
end
