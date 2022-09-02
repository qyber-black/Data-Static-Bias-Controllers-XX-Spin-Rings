function parsave(name, Results, best, fastest, FailedRuns, Info, base, m)
  % Helper function for compute_bias_controls/localisation
  save(['results/data_' name '.mat'], '-v7.3', 'Results', 'best', 'fastest', 'FailedRuns', 'Info');
  n = size(dir(['results/data_' base '-*.mat']),1);
  if m == 0
    if isempty(Results)
      disp(sprintf('Computing %s', name));
    else
      disp(sprintf('Saved data for %s', name));
    end
  else
    if isempty(Results)
      disp(sprintf('Computing %s: %d of %d, %g%% complete', name, n, m, n/m * 100));
    else
      disp(sprintf('Saved data for %s: %d of %d, %g%% complete', name, n, m, n/m * 100));
    end
  end
end
