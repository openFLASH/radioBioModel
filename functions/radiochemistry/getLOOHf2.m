%----------------------------
% The simple heuristic ODE system to describethe effect ofthe pulsed beam
%
% INPUT
% k -_SCALAR VECTOR_- Rate constants of the model
% O2 -_SCALAR_- [oxygen] (mol/l)
% |DRp| -_SCALAR VECTOR_- |DRp(t)| Instantaneous dose rate (Gy/s) at the i-th time point
% |t1| -_SCALAR VECTOR_- Time (ms) at the i-th point
% |no_repainting| -_SCALAR_- Number of time the DRp sequence is repeated, with delay time > 2s
%
% OUTPUT
% |LOOHf| - _SCALAR VECTOR_ - (mol/l) Concentration of LOOH at the end of the experiment
%----------------------------

function LOOHf  = getLOOHf2(k, O2 , DRp , tin , no_repainting, verbose)

      if nargin < 5
        verbose = false;
      end

      TafterPulse = 4; %s Time to continue simulation without beam

      %solve ODE
      y0 = zeros(4,1); %mol/l of the different chemicals at the end of the

      %opts = odeset('NonNegative',1:length(y0),'InitialStep',dt.*1e-4);
      opts = odeset('NonNegative',1:length(y0));

      [t,y] = ode15s(@SimpleModel,[0 , tin(end)] ,y0 , opts  , k, O2 , DRp , tin);


      if verbose
          figs = 1:numel(y0); %List of the species to plot
          [~ , labels] = SimpleModel();
          DRp_t = interp1(tin,DRp,t,'linear',0);
          for fig_i = 1:numel(figs)

            fig = figs(fig_i);
            figure(99+fig)

            taxis = t;
            idx = round(numel(taxis) .* 0.1);
            loglog(taxis,y(:,fig),'.b')
            hold on
            loglog(taxis,DRp_t,'-r') %display beam ON indicator
          end
      end

      %After the last pulse, compute evolution over 2s
      delta = 1e-10; %s When computing decay during beam OFF, stops a bit before next pulse
      Tfinal = t(end) + TafterPulse;
      [t1,y1] = ode15s(@SimpleModel,[t(end)+delta , Tfinal],y(end,:),opts , k, O2 , DRp , tin); %set dose rate =0 to turn beam off

      if verbose
          for fig_i = 1:numel(figs)
            fig = figs(fig_i);
            figure(99+fig)
            hold on
            loglog(t1,y1(:,fig),'-g')
            xlabel('Time (s)')
            ylabel('Concentration (\mu mol/l)')
            title(['[',labels{fig},'] evolution',])
            grid minor
            drawnow
          end
      end

      %Collect the current solution with previous solutions
      t = [t ; t1];
      y = [y ; y1];

      [~ , labels] = SimpleModel();
      LOOHidx = find(strcmp(labels, 'LOOH')); %find the index for [O2]
      LOOHf =  y(end,LOOHidx) .* no_repainting; %M
      %LOOHf =  y(end,LOOHidx); %M
end
