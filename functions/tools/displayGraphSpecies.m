%% displayGraphSpecies
% Display graph for all species
%
%% Syntax
% |displayGraphSpecies(t, y , doseRateIndex , legendSTR3 , colorIndex)|
%
%
%% Description
% |displayGraphSpecies(t, y , doseRateIndex , legendSTR3 , colorIndex)| Description
%
%
%% Input arguments
% |im1| - _STRING_ -  Name
%
%
%% Output arguments
%
% |res| - _STRUCTURE_ -  Description
%
%
%% Contributors

function displayGraphSpecies(t, y , doseRateIndex , legendSTR3 , colorIndex)


    [~ , labels]= radiolysisKinetics2P_a();

    markers = {'-.','-','--',':'};
    colour = {'k','b','g','r','c','m','y'};

    for fig = 1:size(y,2)

      figure(fig)
      %semilogx(t,y(:,fig), markers{doseRateIndex},'MarkerIndices',1:30:length(t),'MarkerSize',5)
      loglog(t,y(:,fig), [markers{mod(doseRateIndex,length(markers))+1} colour{mod(colorIndex,length(colour))+1}])
      hold on
      xlabel('Time (s)')
      ylabel('Concentration (\mu mol/l)')
      title(['[',labels{fig},'] evolution',])
      legend(legendSTR3);
      grid minor

    end



    %pause
    drawnow

end
