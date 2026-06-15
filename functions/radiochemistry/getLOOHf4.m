%----------------------------
% The simple heuristic ODE system to describethe effect ofthe pulsed beam
%
% INPUT
% |k| -_SCALAR VECTOR_- Rate constants of the model
% |O2| -_SCALAR_- [oxygen] (mol/l)
% |totalDose| -_SCALAR_ Delivered dose (Gy)-
% |numberOfPulses| -_SCALAR_ Number of pulses -
% |SweepTime| -_SCALAR_- Time (ms) between 2 pulses
% |DwellTime| -_SCALAR_- Duration time (ms) of 1 pulse
%
% OUTPUT
% |LOOHf| - _SCALAR VECTOR_ - (mol/l) Concentration of LOOH at the end of the experiment
%----------------------------

function LOOHf  = getLOOHf4(k, O2 , totalDose , numberOfPulses , SweepTime , DwellTime, verbose)


      SweepTime = SweepTime .* 1e-3; %Convert to s
      DwellTime = DwellTime .* 1e-3; %Convert to s

      if nargin < 5
        verbose = false;
      end

      TafterPulse = 4; %s Time to continue simulation without beam

      [~ , labels] = SimpleModel();
      LOOHidx = find(strcmp(labels, 'LOOH')); %find the index
      Lridx = find(strcmp(labels, 'L^.'));
      LOOridx = find(strcmp(labels, 'LOO^.'));

      DRp = totalDose ./ (numberOfPulses .* DwellTime); %Instantaneous dose rate

      %solve ODE
      y0 = zeros(4,1); %mol/l of the different chemicals at the end of the

      %opts = odeset('NonNegative',1:length(y0),'InitialStep',dt.*1e-4);
      opts = odeset('NonNegative',1:length(y0));
      t = [];
      y = [];
      t1 = 0;
      ys = y0;

      for pls = 1:numberOfPulses
            %Beam ON
            ts = t1(end);
            te = ts + DwellTime;
            [t1,y1] = ode15s(@SimpleModel, [ts , te] ,ys , opts  , k, O2 , [DRp , DRp] , [ts,te]);
            t = [t;t1];
            y = [y;y1];
            ys = y1(end,:);


            if verbose
                figs = 1:numel(y0); %List of the species to plot
                [~ , labels] = SimpleModel();
                DRp_t = interp1([ts,te],[DRp , DRp],t1,'linear',0);

                for fig_i = 1:numel(figs)
                  fig = figs(fig_i);
                  figure(99+fig)
                  taxis = t1;
                  loglog(taxis(),y1(:,fig),'.b')
                  text(taxis(end),y1(end ,fig),num2str(pls))
                  hold on
                  loglog(taxis,DRp_t,'-r') %display beam ON indicator
                  text(taxis(end),DRp_t(end),num2str(pls))
                end
            end

            %Beam OFF
            ts = t1(end);
            te = ts + SweepTime;
            [t1,y1] = ode15s(@SimpleModel, [ts , te] ,ys , opts  , k, O2 , [0 , 0] , [ts,te]);
            t = [t;t1];
            y = [y;y1];
            ys = y1(end,:);

            if verbose
                figs = 1:numel(y0); %List of the species to plot
                for fig_i = 1:numel(figs)
                  fig = figs(fig_i);
                  figure(99+fig)
                  taxis = t1;
                  loglog(taxis,y1(:,fig),'.g')
                  hold on
                end
            end

      end

      %After the last pulse, compute evolution over 2s
      ts = t(end);
      te = ts + TafterPulse;
      ys = y(end,:);

      [t1,y1] = ode15s(@SimpleModel, [ts , te] ,ys , opts  , k, O2 , [0 , 0] , [ts,te]);
      t = [t;t1];
      y = [y;y1];
      ys = y1(end,:);

      if verbose
          for fig_i = 1:numel(figs)
            fig = figs(fig_i);
            figure(99+fig)
            hold on
            loglog(t1,y1(:,fig),'-r')
            xlabel('Time (s)')
            ylabel('Concentration (\mu mol/l)')
            title(['[',labels{fig},'] evolution',])
            grid minor
            drawnow
          end
      end

      LOOHf =  y(end,LOOHidx); %M
end
