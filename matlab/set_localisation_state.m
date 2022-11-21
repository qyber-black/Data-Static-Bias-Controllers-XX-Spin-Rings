function state = set_localisation_states(id)
  % create_localiation (states)
  %
  % id  Readout id: 1 - .1 readout window
  %                 2 - .1 readout window with .05 pertubation and 500 samples
  %
  % Create static bias control figures for localising spin 1
  % in N rings for N=3:30.

  if ~exist('id','var')
    error('No localisation state specified');
  end

  state.Nmax = 30;           % Max ring size

  state.CT = 1000;           % Localisation time interval (Time = CT/2)
  state.B = NaN;             % Maximum bias?
  state.initT = [];          % (1) - max time, (2) - samples for init times from chain (or empty)
  state.repeats = 2000;      % How many restarts?
  state.symm = 1;            % Enforce symmetry
  state.bias_init = 1;       % Initial biases with peak/trough/constant instead of random (0)?
  state.noise = 0;           % Add noise to initial values (or 0)
  state.min_err = 0.01;      % Largest error acceptable for shortest solution.

  state.readout = [];
  state.id_str = 'unknown';
  state.p_str = 'unknown';

  if id == 1
    state.readout = [state.CT];   % Reatdout window
    state.id_str = 'dt';
    state.p_str = 'T,dT';
  elseif id == 2
    state.readout = [state.CT 0.05 500]; % Readout window (or empty to ignore)
                                   %  readout(1) Time window for readout (0 - do not consider window)
                                   %  readout(2) Pertubation (or 0 for none, time window required)
                                   %  readout(3) Samples for expectation over pertubation (if pertubation not 0)
      state.min_err = 0.1;         % Largest error acceptable for shortest solution.
      state.noise = 200;           % Add noise to initial values (or 0)
      state.id_str = 'dtp';
      state.p_str = 'T,dT,p';
  else
    error(sprintf('Unknown state %g',state));
  end

end