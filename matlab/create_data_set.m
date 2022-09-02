function create_data_set()
  % create_data_set ()
  %
  % Create data set from raw results.

  Nmax = 20;

  for id = [0 1]
    s = set_bias_control_state(id);

    TargetMax = floor((Nmax+1)/2);
    Time = zeros(Nmax,TargetMax);
    Error = zeros(Nmax,TargetMax);
    MeanTime = zeros(1,TargetMax);
    MeanTimeN = zeros(1,TargetMax);
    MaxTime = nan(1,TargetMax);
    MinTime = nan(1,TargetMax);

    for N = 3:Nmax
      for target = 2:floor((N+1)/2) % 1 to target transition (in ring)

        name = sprintf('%s-%d-%d', s.id_str, N, target);

        if exist(['results/data_bias_control_' name '.mat'], 'file')

          % Convert data set for type id, ring size N, 1-target transition
          display(['Converting ' name]);
          data = load(['results/data_bias_control_' name '.mat']);

          info.N = N;
          info.in = data.Info.args.in;
          info.out = data.Info.args.out;
          info.fastest_min_err = data.Info.args.min_err;
          if id == 0
            info.dt = 0;
          else
            info.dt = data.Info.args.readout(1);
          end
          results = data.Results;
          results_idx_best = data.best;
          results_idx_fastest = data.fastest;

          dpdJ = data.Info.args.obj.bias_sensitivity (results,data.Info,'couplings');
          E = arrayfun(@(l) results{l}.err, 1:size(results,2));
          S = arrayfun(@(l) norm(dpdJ{l}), 1:size(results,2));
          [E,IDX] = sort(E);
          S = S(IDX);
          IDX = log10(E) < -1;
          sensitivity.error = E(IDX);
          sensitivity.dpdJ_norm = S(IDX);
          sensitivity.taub = ktaub([E', S'], 0.05);

          save(['data_set/data_' name '.mat'], '-v7.3', 'info', 'results', 'results_idx_best', 'results_idx_fastest', 'sensitivity');

          % Create bias figure
          display(['Bias ' name]);
          figure(1);
          clf;
          plot_bias(info, results, results_idx_best, results_idx_fastest, data.Info.args.obj);
          drawnow();
          refresh();
          savefig(1, ['data_set/bias_' name '.fig'], 'compact');
          screen2png(['data_set/bias_' name]);

          % Create sensitivity figure
          display(['Sensitivity ' name]);
          figure(2);
          clf;
           % Error
          line(1:size(sensitivity.error,2),log10(sensitivity.error),'Color','r','LineWidth',2);
          ax1 = gca;
          ax1.XColor = 'k';
          ax1.YColor = 'r';
          axis tight;
          % Sensitivity
          ax2 = axes('Position',ax1.Position,...
                    'XAxisLocation','top',...
                    'YAxisLocation','right',...
                    'Color','none');
          ax2.XColor = 'k';
          ax2.YColor = 'b';
          line(1:size(sensitivity.dpdJ_norm,2),log10(sensitivity.dpdJ_norm),'Parent',ax2,'Color','b');
          if id == 0
            ylabel(ax2,sprintf('log10(|dp(|out> <- |in>,T)/dJ|)'));
            ylabel(ax1,sprintf('log10(1-p(|out> <- |in>,T))'));
          else
            ylabel(ax2,sprintf('log10(|dp(|out> <- |in>,T)/dJ|)'));
            ylabel(ax1,sprintf('log10(1-p(|out> <- |in>,T))'));
          end
          xlabel(ax1,'Run, ordered by increasing error, arb. units');
          title(sprintf('Error and sensitivity w.r.t. couplings for %d-ring, |1>-|%d>, tau_b=%.3g', N, target, sensitivity.taub));
          axis tight;
          drawnow();
          refresh();
          savefig(2, ['data_set/sensitivity_' name '.fig'], 'compact');
          screen2png(['data_set/sensitivity_' name]);

          if results_idx_fastest > 0
            Time(N,target) = results{results_idx_fastest}.time;
            Error(N,target) = results{results_idx_fastest}.err;
            MeanTime(target) = MeanTime(target) + results{results_idx_fastest}.time;
            MeanTimeN(target) = MeanTimeN(target) + 1;
            MinTime(target) = min(MinTime(target),results{results_idx_fastest}.time);
            MaxTime(target) = max(MaxTime(target),results{results_idx_fastest}.time);
          end

        else
          error(['Missing ' name]);
        end

      end
    end

    % Summary plot for fastest results
    display(['Fastest ' s.id_str]);
    figure(3);
    clf;
    h = bar3(Time(3:Nmax,2:TargetMax));
    for col = 1:TargetMax-1
      zdata = [];
      for row = 1:Nmax-2
        zdata = [zdata; ones(6,4)*Error(row,col)];
      end
      h(col).CData = zdata;
    end
    remove_empty_bars(h);
    h = colorbar;
    set(get(h,'ylabel'),'String', 'Infidelity');
    axis tight;
    xlabel('Target state |n>');
    ax = gca;
    ax.XTick = 1:TargetMax-1;
    ax.XTickLabel = 2:TargetMax;
    ylabel('Ring size N');
    ax.YTick = 1:Nmax-2;
    ax.YTickLabel = 3:Nmax;
    zlabel('Fastest time');
    title('Fastest times for state transfer from |1> to |n> for N-rings');
    hold on;
    h = plot3(1:TargetMax-1,zeros(1,TargetMax-1),MeanTime(2:TargetMax)./MeanTimeN(2:TargetMax),'.r');
    set(h, 'MarkerSize', 25);
    for l = 2:TargetMax
      h=plot3([l-1; l-1], [0; 0], [MinTime(l); MaxTime(l)], '-r');
      set(h, 'LineWidth', 2);
    end
    drawnow ();
    refresh ();
    % Save
    savefig(3, ['data_set/fastest_' s.id_str '.fig'], 'compact');
    screen2png(['data_set/fastest_' s.id_str]);

  end

  % Localisation
  s = set_localisation_state(1);
  for N = 3:Nmax       % Ring size
    name = sprintf('dt-%d-1', N);
    fname = sprintf('results/data_localisation_dt-%d.mat', N);
    % Find bias controls
    if exist(fname, 'file')

      % Converting localisation data
      display(['Converting ' name]);
      data = load(fname);

      info.N = N;
      info.in = data.Info.args.in;
      info.out = data.Info.args.out;
      info.dt = data.Info.args.readout(1);
      results = data.Results;
      results_idx_best = data.best;

      dpdJ = data.Info.args.obj.bias_sensitivity (results,data.Info,'couplings');
      E = arrayfun(@(l) results{l}.err, 1:size(results,2));
      S = arrayfun(@(l) norm(dpdJ{l}), 1:size(results,2));
      [E,IDX] = sort(E);
      S = S(IDX);
      IDX = log10(E) < -1;
      sensitivity.error = E(IDX);
      sensitivity.dpdJ_norm = S(IDX);
      % Kendall tau_b
      sensitivity.taub = ktaub([E', S'], 0.05);

      save(['data_set/data_' name '.mat'], '-v7.3', 'info', 'results', 'results_idx_best', 'sensitivity');

      % Create bias figure
      display(['Bias ' name]);
      figure(1);
      clf;
      tt = info.dt;
      plot_time = [0:tt/1000:tt];
      % natural evolution
      plot_nat = data.Info.args.obj.trace(info.in,info.out,plot_time);
      % evolution with control
      H = data.Info.args.obj.H + diag(results{results_idx_best}.bias);
      [V,e] = eig (H);
      e = diag(e);
      E = cellfun (@(x) V * diag(exp(-i * x * e)) * V', num2cell (plot_time), 'UniformOutput', false);
      plot_ctrl =  cellfun (@(x) abs(x(info.out, info.in))^2, E);
      plot_result(results_idx_best, plot_time, plot_ctrl, plot_nat, 3,2,1,3,'Best solution');
      subplot(3,2,1);
      title(sprintf('Bias control for |1>-|1> transition in ring of N=%d spins', info.N));
      subplot(3,2,[5 6]);
      plot_eigenstructure([V;e'],'best solution');
      % Error Histogram
      subplot(3,2,[2 4]);
      Err = log10(arrayfun(@(x) (results{x}.err),[1:size(results,2)]));
      histogram(Err,'FaceColor',[0 0 1],'FaceAlpha',1);
      title(sprintf('log(Error) Histogram over %d runs', size(results,2)));
      axis tight;
      drawnow();
      refresh();
      savefig(1, ['data_set/bias_' name '.fig'], 'compact');
      screen2png(['data_set/bias_' name]);

      % Create sensitivity figure
      display(['Sensitivity ' name]);
      figure(2);
      clf;
      % Error
      line(1:size(sensitivity.error,2),log10(sensitivity.error),'Color','r','LineWidth',2);
      ax1 = gca;
      ax1.XColor = 'k';
      ax1.YColor = 'r';
      axis tight;
      % Sensitivity
      ax2 = axes('Position',ax1.Position,...
                'XAxisLocation','top',...
                'YAxisLocation','right',...
                'Color','none');
      ax2.XColor = 'k';
      ax2.YColor = 'b';
      line(1:size(sensitivity.dpdJ_norm,2),log10(sensitivity.dpdJ_norm),'Parent',ax2,'Color','b');
      ylabel(ax2,sprintf('log10(|dp(|out> <- |in>,T)/dJ|)'));
      ylabel(ax1,sprintf('log10(1-p(|out> <- |in>,T))'));
      xlabel(ax1,'Run, ordered by increasing error, arb. units');
      title(sprintf('Error and sensitivity w.r.t. couplings for %d-ring, |1>-|1>, tau_b=%.3g', N, sensitivity.taub));
      axis tight;
      drawnow();
      refresh();
      savefig(2, ['data_set/sensitivity_' name '.fig'], 'compact');
      screen2png(['data_set/sensitivity_' name]);

    else
      error(['Missing ' name]);
    end
  end

  function plot_bias (info, results, results_idx_best, results_idx_fastest, qsn)
    if results_idx_best > 0
      % Setup trace plots
      tt = results{results_idx_best}.time + info.dt/2;
      plot_time = [0:tt/1000:tt];
      % natural evolution
      plot_nat = qsn.trace(info.in,info.out,plot_time);
      % evolution with control
      H = qsn.H + diag(results{results_idx_best}.bias);
      [V,e] = eig (H);
      e = diag(e);
      E = cellfun (@(x) V * diag(exp(-i * x * e)) * V', num2cell (plot_time), 'UniformOutput', false);
      plot_ctrl =  cellfun (@(x) abs(x(info.out, info.in))^2, E);
      plot_result(results_idx_best, plot_time, plot_ctrl, plot_nat, 3,3,1,4,'Best solution');
      subplot(3,2,5);
      plot_eigenstructure([V;e'],'best solution');
    end
    % Fastest solution
    if results_idx_fastest > 0
      % Setup trace plots
      tt = results{results_idx_fastest}.time + info.dt/2;
      plot_time = [0:tt/1000:tt];
      % natural evolution
      plot_nat = qsn.trace(info.in,info.out,plot_time);
      % evolution with control
      H = qsn.H + diag(results{results_idx_fastest}.bias);
      [V,e] = eig (H);
      e = diag(e);
      E = cellfun (@(x) V * diag(exp(-i * x * e)) * V', num2cell (plot_time), 'UniformOutput', false);
      plot_ctrl =  cellfun (@(x) abs(x(info.out, info.in))^2, E);
      plot_result(results_idx_fastest, plot_time, plot_ctrl, plot_nat, 3,3,2,5,'Fastest solution');
      subplot(3,2,6);
      plot_eigenstructure([V;e'],'fastest solution');
    end
    % Error vs time
    subplot(3,3,3);
    times = arrayfun(@(x) (results{x}.time),[1:size(results,2)]);
    plot(log10(arrayfun(@(x) (results{x}.err),[1:size(results,2)])), times, '*b');
    xlabel('log(Error)');
    ylabel('Time');
    ylim([max(0,min(times)) max(times)]);
    axis tight;
    % Error Histogram
    figure(1);
    subplot(3,3,6);
    Err = log10(arrayfun(@(x) (results{x}.err),[1:size(results,2)]));
    histogram(Err,'FaceColor',[0 0 1],'FaceAlpha',1);
    title(sprintf('Log(Error) histogram over %d runs', size(results,2)));
    axis tight;
    mtit(sprintf('Bias control for |%d>-|%d> transition in ring of N=%d spins', info.in, info.out, info.N));
  end

  % Plot results
  function plot_result (run, plot_time, plot_ctrl, plot_nat, X, Y, trace_fig, bias_fig, str)
    % Plot traces
    subplot(X,Y,trace_fig);
    plot(plot_time, plot_ctrl, plot_time, plot_nat);
    axis([0 plot_time(end) 0 1]);
    %legend('control','no control', 'Location', 'best')
    xlabel('Time t in 1/J')
    ylabel(sprintf('Probability p(|%i> <- |%i>,t)',info.out,info.in))
    axis tight;
    % Plot bias
    subplot(X,Y,bias_fig);
    bar(1:info.N,results{run}.bias,'b');
    xlabel('Spin #');
    ylabel('Bias');
    title(sprintf('%s, T=%.6g, Error=%.6g',str, results{run}.time, results{run}.err));
    axis tight;
  end

  function plot_eigenstructure(V,str)
    % Eigenstructure plot
    colormap(flipud(gray));
    VV = abs(V);
    clim = [min(min(VV(1:end-1,1:end))) max(max(VV(1:end-1,1:end)))];
    im = imagesc(VV,clim);
    imAxes =get(im,'parent');
    hFig = get(imAxes,'parent');
    fs = getBestFontSize(imAxes);
    axis off;
    [rows, cols] = size(V);
    midValue = mean(get(gca,'CLim'));
    ci = (VV < midValue) + 1;
    cmap = colormap ();
    [mx my] = size(cmap);
    cmap = [cmap(1,:); cmap(mx,:)];
    textHandles = zeros(size(V))';
    for i = 1:rows
      for j = 1:cols
        c = cmap(ci(i,j),:);
        if i == info.in
          c = [0 ((2*c(1)-1)/4+3/4) 0];
        elseif i == info.out
          c = [((2*c(1)-1)/4+3/4) 0 0];
        elseif i == rows
          c = ((2*c(1)-1)/4+3/4) * [1 0 1];
        else
          c = ((2*c(1)-1)/4+3/4) * [0 1 1];
        end
        textHandles(j,i) = text(j,i,num2str(V(i,j),'%.3g'),...
                                'color', c,...
                                'horizontalAlignment','center','verticalAlignment','middle',...
                                'fontsize',fs,'clipping','on','visible','on');
      end
    end
    set(imAxes,'UserData',textHandles);
    title(['Eigenstructure of ' str]);
  end

  function fs = getBestFontSize(imAxes)
    % Try to keep font size reasonable for text
    hFig = get(imAxes,'Parent');
    magicNumber = 80;
    nrows = diff(get(imAxes,'YLim'));
    ncols = diff(get(imAxes,'XLim'));
    if ncols < magicNumber && nrows < magicNumber
      ratio = max(get(hFig,'Position').*[0 0 0 1])/max(nrows,ncols);
    elseif ncols < magicNumber
      ratio = max(get(hFig,'Position').*[0 0 0 1])/ncols;
    elseif nrows < magicNumber
      ratio = max(get(hFig,'Position').*[0 0 0 1])/nrows;
    else
      ratio = 1;
    end
    fs = min(10,ceil(ratio/5));    % the gold formula
  end

  function remove_empty_bars(hBars)
    % Remove 0 bars from bar3 plot
    for iSeries = 1:numel(hBars)
      zData = get(hBars(iSeries),'ZData');  % Get the z data
      index = logical(kron(zData(2:6:end,2) == 0,ones(6,1)));  % Find empty bars
      zData(index,:) = nan;                 % Set the z data for empty bars to nan
      set(hBars(iSeries),'ZData',zData);    % Update the graphics objects
    end
  end

  function screen2png(filename)
    % Save current figure into filename as png
    oldscreenunits = get(gcf,'Units');
    oldpaperunits = get(gcf,'PaperUnits');
    oldpaperpos = get(gcf,'PaperPosition');
    oldouterpos = get(gcf,'outerposition');
    set(gcf,'units','normalized','outerposition',[0 0 1 1]);
    set(gcf,'Units','pixels');
    scrpos = get(gcf,'Position');
    newpos = scrpos/100;
    set(gcf,'PaperUnits','inches', 'PaperPosition',newpos)
    print('-dpng', filename, '-r100');
    drawnow
    set(gcf,'Units',oldscreenunits, 'PaperUnits',oldpaperunits, 'PaperPosition',oldpaperpos, 'outerposition', oldouterpos);
    drawnow
  end
end
