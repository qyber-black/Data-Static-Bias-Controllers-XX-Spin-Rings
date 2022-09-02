function compute_localisation(states)
  % compute_localiation (states)
  %
  % states  Readout id: 1 - .1 readout window
  %                     2 - .1 readout window with .05 pertubation and 500 samples
  %         Can be a vector for multiple options
  %
  % Create static bias control figures for localising spin 1
  % in N rings for N=3:30.

  if ~exist('states','var')
    states = [1 2];
  end

  problems = [];
  for id = states
    s = set_localisation_state(id);
    for N = 3:s.Nmax % Ring size
      problems = [problems; id, N];
    end
  end
  M = size(problems,1);
  parfor l = 1:M
    id = problems(l,1);
    N = problems(l,2);
    % Setup network
    s = set_localisation_state (id);
    ring = qsn.QSN ('ring', N);
    name = sprintf('localisation_%s-%d', s.id_str, N);
    % Find bias controls
    if ~exist(['results/data_' name '.mat'], 'file')
      parsave(name, [], [], [],[], [], sprintf('localisation_%s',s.id_str) , size(find(problems(:,1) == id),1));
      [Results,best,fastest,FailedRuns,info] = ring.bias_control(1,1, s.CT/2, s.B, s.readout, s.min_err, s.repeats, s.symm, s.initT, s.bias_init, s.noise, 0);
      parsave(name, Results, best, fastest, FailedRuns, info, sprintf('localisation_%s',s.id_str) , size(find(problems(:,1) == id),1));
    else
      disp(sprintf('%s: OK', name));
    end
  end

end
