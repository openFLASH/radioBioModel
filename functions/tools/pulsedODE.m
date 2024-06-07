%% pulsedODE
% Solve the kinetic system for a pulsed beam by sequentially calling ODE
% 1) during the beam on
% 2) during the beam off
% 3) After the last beam off, to describe the final decay
%
% This function avoids using MaxStep to make sure the solver does not step over some behavior that occurs only once in the integration interval.
% The interval is broken into two pieces (during pulse and during beam off time in between pulses) and the solver is called twice.
%
%
%% Syntax
% |[t,y] = pulsedODE(odefcnComprehensive,t0,y0,opts,param)|
%
%
%% Description
% |[t,y] = pulsedODE(odefcnComprehensive,t0,y0,opts,param)| Description
%
%
%% Input arguments
% |@radioModel| -_FUNCTION POINTER_- POinter to the radio kinetic model function
%
% |y0| - _SCALAR VECTOR_ - Inital concentration of the chemica lspecies of the model at 1us post irradiation
%
% |t0| -_SCALAR_- Time (s) at which the |y0| are computed
%
% |opt| -_STRUCT_- Optinal parameter for ode functions
%
% |param| - _STRUCTURE_ - List of parameters for |@radioModel|
%
%% Output arguments
%
% |res| - _STRUCTURE_ -  Description
%
%
%% Contributors
% Authors : Rudi Labarbe (open.reggui@gmail.com)

function [t,y] = pulsedODE(varargin)

  odefcnComprehensive = varargin{1};
  t0 = varargin{2};
  y0 = varargin{3};
  opts = varargin{4};
  param = varargin{5};

  %param.R0   : Average dose rate
  %param.T    : s %Period of the pulsed beam
  %param.t_on : s Beam on time during one period = pulse duration
  %param.td =  s Total Beam ON time to deliver the full dose
  % TimeSimuMaxLocal = s Total duration of the simulation

  figs = [1,8,9,10]; %List of the species to plot

  NbPulses = ceil(param.td ./ param.T); %Number of pulses
  Tstart = t0; %the first pulse starts at t=0s
  Tend = param.t_on;
  Tstart0 = 0;
  t = [];
  y = [];
  delta = 1e-10; %s When computing decay during beam OFF, stops a bit before next pulse

  if param.t_on < 1e-6
    %The model describes what happens after the begining of the homogeneous chemical phase
    % This is after 1us
    %If the BEAM ON time is so short, the model is probably not appropriate
    error('Beam ON time less than 1us')
  end

  %Loop for every pulse
  for puls = 1:NbPulses

    %Solve ODE during BEAM ON
    [t1,y1] = ode15s(odefcnComprehensive,[Tstart , Tend],y0,opts,varargin{5:end});

    %Collect the current solution with previous solutions
    t = [t ; t1];
    y = [y ; y1];


% [~ , labels]= radiolysisKinetics2P_a();
% for fig_i = 1:numel(figs)
%
%   ymin = 1e-4;
%   fig = figs(fig_i);
%   figure(99+fig)
%
%   taxis =t1-Tstart0;
%   idx = round(numel(taxis) .* 0.1);
%   loglog(taxis,y1(:,fig),'.b')
%   text(taxis(idx),y1(idx,fig),num2str(puls))
%   hold on
%   beam = IsBeamOn(t, param) .* ymin .* 1.5;
%   taxis =t-Tstart0;
%   loglog([1e-10,taxis(end)],[beam(1),beam(end)],'-r') %display beam ON indicator
%   xlabel('Time (s)')
%   ylabel('Concentration (\mu mol/l)')
%   title(['[',labels{fig},'] evolution',])
%
% end


    %Run computation during BEAM OFF
    Tstart = t1(end);
    Tend = puls .* param.T - delta; %Stop computation a litle bit before the start of next pulse
    y0 = y(end,:); %Start from the last place where we stopped

    [t1,y1] = ode15s(odefcnComprehensive,[Tstart , Tend],y0,opts,varargin{5:end});

    %Collect the current solution with previous solutions
    %The first time point is the same as last time point of previous computation
    %The last time point will be the same as the first time point of next computation
    t = [t ; t1];
    y = [y ; y1];

    %Preapare next round. Define the timing from start to end of beam ON
    %Tstart = puls .* param.T;
    Tstart = t1(end);
    Tend = Tstart + param.t_on;
    y0 = y(end,:); %Start from the last place where we stopped

% for fig_i = 1:numel(figs)
%   fig = figs(fig_i);
%   figure(99+fig)
%   loglog(t1-Tstart0,y1(:,fig),'-g')
%   hold on
%   xlabel('Time (s)')
%   ylabel('Concentration (\mu mol/l)')
%   title(['[',labels{fig},'] evolution',])
%   newLim = ylim();
%   newLim(1) = ymin;
%   ylim(newLim)
%   grid minor
%   drawnow
% end
% Tstart0 = Tstart;

  end

%After the last pulse, compute evolution over 2s
Tend = Tstart + 2;
[t1,y1] = ode15s(odefcnComprehensive,[Tstart , Tend],y0,opts,varargin{5:end});

%Collect the current solution with previous solutions
t = [t ; t1];
y = [y ; y1];


%pause

end
