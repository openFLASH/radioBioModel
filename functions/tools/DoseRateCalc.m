%% [DRp , t] = DoseRateCalc(no_spots_X , distance_X , switch_X , beam_on_spot)
%
% Compute the instantaneous dose rate as a function of time
%
%% INPUT
% |totalDose| -_SCALAR_- Total delivered dose (Gy)
% |no_spots| -_SCALAR_- Total number of spot (Nb spots in spot map * no_repainting)
% |irradiatedArea| -_SCALAR_- Irradiated field area (mm2)
% |switch_X| -_SCALAR_- Sweep time (ms) between neighbouring spot (along a line or column, assumed identical)
% |beam_on_spot| -_SCALAR_- Dwell time (ms) on one spot (assumed same for all spots)
% |btw_field_time| -_SCALAR_- Time (ms) beswitch_Xtween 2 repainting
% |no_repainting| -_SCALAR_- Number of repainting
%
%% OUTPUT
% |DRp| -_SCALAR VECTOR_- |DRp(t)| Instantaneous dose rate (Gy/s) at the i-th time point at isocenter
% |t| -_SCALAR VECTOR_- Time (s) at the i-th point
% |cumulative_dose_at_point| -_SCALAR VECTOR_- cumulative dose (Gy) at i-th time point
% |pbs_dr_spot| -_SCALAR_- PBS dose rate at isocenter
%
function [DRp , t , cumulative_dose_at_point , pbs_dr_spot] = DoseRateCalc(totalDose , no_spots , irradiatedArea , switch_X , beam_on_spot , btw_field_time , no_repainting )

  if btw_field_time < 2000
    %time between repainting is less than 2s
    %Compute one long string of instantaneous dose rate
    [total_time_for_spot, cumulative_dose_at_point, pbs_dr_spot] = getCumDoseAtISO(no_spots , irradiatedArea , switch_X , beam_on_spot, btw_field_time , no_repainting);
  else
    %Time between repainting is more than 2s
    %compute one single instantaneous DR and then repeat it by the number of repainting
    [total_time_for_spot, cumulative_dose_at_point, pbs_dr_spot] = getCumDoseAtISO(no_spots./no_repainting , irradiatedArea , switch_X , beam_on_spot, 0 , 1);
    totalDose = totalDose ./ no_repainting;

  end

  %Rescale for the total dose
  pbs_dr_spot = pbs_dr_spot .* totalDose ./ cumulative_dose_at_point(end);
  cumulative_dose_at_point =  cumulative_dose_at_point .* totalDose ./ cumulative_dose_at_point(end);

  DRp = diff(cumulative_dose_at_point) ./ (diff(total_time_for_spot).*1e-3); %Instantaneous dose rate (Gy/s)
  t = total_time_for_spot(2:end).*1e-3;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculate PBS dose rate
% 0) Decide what to plot and calculate
% 1) give the required parameters (beam model, grid size, timing info ...)
% and choose the point of interest for the PBS dose rate
% 2) calculate the dose distribution
% 3) calculate the PBS dose rate (ADR)
% 4) calculate the Van de Water dose rate (DADR)
% 5) calculate the field dose rate (dose/irradiation time)

%REFERENCE:
% 1) Toschini, M., Psoroulas, S., Colizzi, I. & Lomax, A. J. Medical physics dataset article : A database of FLASH murine in vivo studies. Med. Phys. 52, 5115–5123 (2025).
% 2) Colizzi, I., Toschini, M., Lomax, A. J. & Psoroulas, S. Systematic analysis of biological endpoint variability and implications for quantitative modeling of the FLASH sparing effect. Phys. Imaging Radiat. Oncol. 37, 100915 (2026).
% 3) https://zenodo.org/records/10886631
% 4) https://exuberant-beak-513.notion.site/354aa5e33bad4771b65168bb613c5768?v=d1046762681f4605af9f4ab756532030&pvs=4
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [total_time_for_spot, cumulative_dose_at_point , pbs_dr_spot] = getCumDoseAtISO(no_spots_in , irradiatedArea , switch_Xin , beam_on_spotin , btw_field_time , no_repainting)
      %% 0) decide what to plot and calculate

      dose_distribtuion = false; %plot dose distribution

      vdw_dose_rate = false; % mean DADR over all voxels
      vdw_dose_rate_spot_in_target = false; % mean DADR over all spots
      field_dose_rate = false; % field dose rate

      pbs_dose_rate = false; % mean ADR over all voxels
      pbs_dose_rate_spot_in_target = false; % mean ADR over all spots
      pbs_dose_rate_point_of_interest = true; % ADR at point of interest


      %% 1) define the parameters
      % define gaussian beam model
      sigma_X = 0.58; %cm
      sigma_Y = 0.58; %cm
      weight = 1;

      % define weights for spots at the rim and corners (default = 1)
      rim_weight = 1;
      corner_weight = 1;

      % define grid
      [~ , no_spots_X , no_spots_Y] = multipliers(no_spots_in./no_repainting); %divide the number of spot to form the rectangle closest to a square

      distance_X = 0.1 .* sqrt(irradiatedArea ./ ((no_spots_X-1) .* (no_spots_Y-1))); %cm Distance between spots
      distance_Y = 0.1 .* sqrt(irradiatedArea ./ ((no_spots_X-1) .* (no_spots_Y-1)));

      % define spot changing time
      switch_X = switch_Xin;
      switch_Y = switch_Xin;
      beam_on_spot = beam_on_spotin;

      % point of interest for the PBS dose rate
      point_of_interest = [distance_X.*(no_spots_X-1)./2, distance_Y.*(no_spots_Y-1)./2]; %Center of irradiated field

      % dose level for the irradiated area
      level = 0.95;

      % dose-threshold for folkers
      threshold = 0.05;

      % spots inside the target
      %ordered_points_in_target =  [ 7 8 9 12 13 14 17 18 19 22 23 24 ];

      % repeanting
      % btw_field_time = 0;
      % no_repainting = 1;

      %% 2) order spots and plot dose distribution

      % create grid
      positions_X = linspace(0,(no_spots_X -1)*distance_X,no_spots_X);
      positions_Y = linspace(0,(no_spots_Y-1)*distance_Y,no_spots_Y);
      [x,y] = meshgrid(positions_X,positions_Y);

      % create fine grid for plotting
      [xx, yy] = meshgrid(positions_X(1)-1.5:0.05:positions_X(end)+1.5, positions_Y(1)-1.5:0.05:positions_Y(end)+1.5);

      % order the spot in an "S" pattern
      ordered_points = zeros(numel(x), 2);
      index = 1;
      for row = 1:size(x, 1)
          % Alternate the direction along X for each row
          if mod(row, 2) == 1
              % Traverse from left to right along X
              ordered_points(index:index + size(x, 2) - 1, :) = [x(row, :)', y(row, :)'];
          else
              % Traverse from right to left along X
              ordered_points(index:index + size(x, 2) - 1, :) = [x(row, end:-1:1)', y(row, end:-1:1)'];
          end
          index = index + size(x, 2);
      end
      ordered_points = flipud(ordered_points);

      %% Compute Dose Distribution
      dose_distribution = zeros(size(xx));
      for i = 1:size(ordered_points, 1)
          position = ordered_points(i, :);
          if (position(1) == 0 || position(1) == (no_spots_X-1)  * distance_X) && (position(2) == 0 || position(2) == (no_spots_Y-1)  * distance_Y)
              weight_factor = corner_weight;
          elseif position(1) == 0 || position(1) == (no_spots_X-1)  * distance_X || position(2) == 0 || position(2) == (no_spots_Y-1)  * distance_Y
              weight_factor = rim_weight;
          else
              weight_factor = 1;
          end
          dose_distribution = dose_distribution + weight_factor * weight * ...
              exp(-((xx - ordered_points(i, 1)).^2 / (2 * sigma_X^2) + ...
              (yy - ordered_points(i, 2)).^2 / (2 * sigma_Y^2)));
      end

      %% Define the Region of Interest (ROI)
      mask = (xx >= 0 & xx <= (no_spots_X-1) * distance_X) & (yy >= 0 & yy <= (no_spots_Y-1) * distance_Y);
      dose_distribution(~mask) = NaN; % Set outside region to NaN

      %% Plot Dose Distribution in Defined Region
      if dose_distribtuion == true
          figure;
          contourf(xx, yy, dose_distribution, 200, 'LineColor', 'none');
          title('Dose Distribution in Defined Region');
          grid on
          xlabel('X-axis (cm)');
          ylabel('Y-axis (cm)');
          colorbar;
          colormap('jet');
          hold on;

          % Plot spot positions
          plot(ordered_points(:,1), ordered_points(:,2), ':k', 'MarkerEdgeColor', 'k');
          scatter(ordered_points(:,1), ordered_points(:,2), 'k', 'MarkerEdgeColor', 'k');
          scatter(ordered_points(1,1), ordered_points(1,2), 'rd', 'filled', 'MarkerEdgeColor', 'k');
          scatter(ordered_points(end,1), ordered_points(end,2), 'rx', 'filled', 'MarkerEdgeColor', 'k');

          % Plot point of interest
          scatter(point_of_interest(1), point_of_interest(2), '*', 'filled', 'MarkerEdgeColor', 'b');

          % Plot dose contour
          contour(xx, yy, dose_distribution, [level*max(dose_distribution(:)), level*max(dose_distribution(:))], 'LineColor', 'k');

          hold off;
      end

      %% Compute Dose Rate in the Defined Region
      time = (switch_X*(no_spots_X-1)  + switch_Y*(no_spots_Y-1)  + ...
          max(dose_distribution(:))/no_repainting * beam_on_spot * (no_spots_Y+no_spots_X-2)) ...
          * no_repainting + btw_field_time * (no_repainting - 1);

      dose_rate_distribution = dose_distribution ./ (0.001 * time);
      dose_rate_distribution(~mask) = NaN; % Mask out unwanted regions

      %% Plot Dose Rate Distribution (Limited to Defined Region)
      if(0)
          figure;
          contourf(xx, yy, dose_rate_distribution, 100, 'LineColor', 'none');
          colormap('jet');
          colorbar;
          hold on;

          % Plot spot positions
          plot(ordered_points(:,1), ordered_points(:,2), ':k', 'MarkerEdgeColor', 'k');
          scatter(ordered_points(:,1), ordered_points(:,2), 'k', 'MarkerEdgeColor', 'k');
          scatter(ordered_points(end,1), ordered_points(end,2), 'rx', 'filled', 'MarkerEdgeColor', 'k');

          % Plot point of interest
          scatter(point_of_interest(1), point_of_interest(2), '*', 'filled', 'MarkerEdgeColor', 'b');

          % Plot dose contour
          contour(xx, yy, dose_rate_distribution, [level*max(dose_rate_distribution(:)), level*max(dose_rate_distribution(:))], 'LineColor', 'k');

          hold off;
          grid on;
          xlabel('X Position (cm)');
          ylabel('Y Position (cm)');
          title('Field Dose Rate Distribution (Gy/s)');
          set(gca, 'YDir', 'normal');
        end

      %% 3a) calcualte PBS dose rate at the point of interest

      if pbs_dose_rate_point_of_interest == true

          cumulative_dose_at_point = zeros(size(ordered_points, 1)*no_repainting, 1);
          total_time_for_spot = zeros(size(ordered_points, 1)*no_repainting, 1);

          for repeanting = 1:no_repainting

              if repeanting == 1
                  total_time_for_spot(1) = beam_on_spot * corner_weight * weight/no_repainting;
                  cumulative_dose_at_point(1) =  corner_weight * weight/no_repainting * exp(-((point_of_interest(1) - ordered_points(1, 1)).^2 / (2 * sigma_X^2) + (point_of_interest(2) - ordered_points(1, 2)).^2 / (2 * sigma_Y^2)));
              end

              for i = 2:size(ordered_points, 1)

                  % calculate cumulative dose
                  position = ordered_points(i, :);
                  if (position(1) == 0 || position(1) == (no_spots_X-1)  * distance_X) && (position(2) == 0 || position(2) == (no_spots_Y-1)  * distance_Y)
                      weight_factor = corner_weight;
                  elseif position(1) == 0 || position(1) == (no_spots_X-1)  * distance_X || position(2) == 0 || position(2) == (no_spots_Y-1)  * distance_Y
                      weight_factor = rim_weight;
                  else
                      weight_factor = 1;
                  end

                  dose_at_point_of_interest = weight_factor * weight/no_repainting * exp(-((point_of_interest(1) - ordered_points(i, 1)).^2 / (2 * sigma_X^2) + (point_of_interest(2) - ordered_points(i, 2)).^2 / (2 * sigma_Y^2)));
                  cumulative_dose_at_point(i+ size(ordered_points, 1)*(repeanting-1)) =  cumulative_dose_at_point(i-1 + size(ordered_points, 1)*(repeanting-1)) + dose_at_point_of_interest;
                  switch_time_T = switch_X * (mod(i-1, (no_spots_X) ) > 0);
                  switch_time_U = switch_Y * (mod(i-1, (no_spots_Y) ) == 0 && i > 1);
                  total_time_for_spot(i+ size(ordered_points, 1)*(repeanting-1)) = total_time_for_spot(i-1 + size(ordered_points, 1)*(repeanting-1)) + switch_time_T + switch_time_U + weight_factor * weight/no_repainting * beam_on_spot;

              end

              if no_repainting > 1 & no_repainting > repeanting
                  total_time_for_spot(size(ordered_points, 1)*repeanting + 1 ) = total_time_for_spot(size(ordered_points, 1)*repeanting ) + btw_field_time;
                  cumulative_dose_at_point(size(ordered_points, 1)*repeanting + 1 ) = cumulative_dose_at_point(size(ordered_points, 1)*repeanting );
              end

          end

          cumulative_dose_at_point = [0.0; cumulative_dose_at_point];
          total_time_for_spot = [ 0.0; total_time_for_spot];

          % Find the treshold
          upper_threshold = 1 - threshold/2;
          lower_threshold = threshold/2;

          % find the time within the treshold of the total dose
          target_dose_down_percent = lower_threshold * cumulative_dose_at_point(end);
          target_dose_up_percent = upper_threshold * cumulative_dose_at_point(end);
          cumulative_dose_at_point = round(cumulative_dose_at_point, 4);
          [unique_dose, unique_indices] = unique(cumulative_dose_at_point);
          unique_total_time = total_time_for_spot(unique_indices);
          time_at_down_percent = interp1(unique_dose, unique_total_time, target_dose_down_percent);
          time_at_up_percent = interp1(unique_dose, unique_total_time, target_dose_up_percent);
          dose_at_down_percent = interp1(total_time_for_spot, cumulative_dose_at_point, time_at_down_percent);
          dose_at_up_percent = interp1(total_time_for_spot, cumulative_dose_at_point, time_at_up_percent);

          time_interval = time_at_up_percent - time_at_down_percent;
          dose_interval = dose_at_up_percent - dose_at_down_percent;

          pbs_dr_spot = dose_interval/(0.001*time_interval);
          if(0)
              disp(' ');
              disp(['---- PBS dose rate -----']);
              disp(['dose threshold: ', num2str(threshold*100), '%']);
              disp(['time within interval: ', num2str(time_interval), ' ms']);
              disp(['dose within interval: ', num2str(dose_interval), ' Gy']);

              disp(['pbs dose rate: ', num2str(pbs_dr_spot), ' Gy/s']);

              % Plot the cumulative dose at the point of interest as a function of time
              figure;
              plot(total_time_for_spot, cumulative_dose_at_point, 'o-');
              title('Cumulative Dose at a Point as a Function of Time');
              xlabel('Time (ms)');
              ylabel('Cumulative Dose at Point of Interest (Gy)');
              grid on;

          end
      end

      %% 3b) calcualte PBS dose rate at all point of interest

      if pbs_dose_rate_spot_in_target == true

          point_of_interest = ordered_points;

          figure;
          disp(' ');
          disp('---- PBS dose rate for each spot-----');
          hold on

          for point = 2:size(point_of_interest,1)

              cumulative_dose_at_point = zeros(size(ordered_points, 1)*no_repainting, 1);
              total_time_for_spot = zeros(size(ordered_points, 1)*no_repainting, 1);

              for repeanting = 1:no_repainting

                  % first spot
                  if repeanting == 1
                      total_time_for_spot(1) = beam_on_spot * corner_weight * weight/no_repainting;
                      cumulative_dose_at_point(1) =  corner_weight * weight/no_repainting * exp(-((point_of_interest(1) - ordered_points(1, 1)).^2 / (2 * sigma_X^2) + (point_of_interest(2) - ordered_points(1, 2)).^2 / (2 * sigma_Y^2)));
                  end

                  for i = 2:size(ordered_points, 1)

                      position = ordered_points(i-1, :);
                      if (position(1) == 0 || position(1) == (no_spots_X-1)  * distance_X) && (position(2) == 0 || position(2) == (no_spots_Y-1)  * distance_Y)
                          weight_factor = corner_weight;
                      elseif position(1) == 0 || position(1) == (no_spots_X-1)  * distance_X || position(2) == 0 || position(2) == (no_spots_Y-1)  * distance_Y
                          weight_factor = rim_weight;
                      else
                          weight_factor = 1;
                      end

                      dose_at_point_of_interest = weight_factor * weight/no_repainting * exp(-((point_of_interest(point,1) - ordered_points(i, 1)).^2 / (2 * sigma_X^2) + (point_of_interest(point,2) - ordered_points(i, 2)).^2 / (2 * sigma_Y^2)));
                      cumulative_dose_at_point(i + size(ordered_points, 1)*(repeanting-1)) = cumulative_dose_at_point(i-1 + size(ordered_points, 1)*(repeanting-1)) + dose_at_point_of_interest;

                      switch_time_T = switch_X * (mod(i-1, (no_spots_X) ) > 0);
                      switch_time_U = switch_Y * (mod(i-1, (no_spots_Y) ) == 0 && i > 1);
                      total_time_for_spot(i + size(ordered_points, 1)*(repeanting-1)) = total_time_for_spot(i-1 + size(ordered_points, 1)*(repeanting-1)) + switch_time_T + switch_time_U + weight_factor * weight/no_repainting * beam_on_spot;

                  end

                  if no_repainting > 1 & no_repainting > repeanting
                      total_time_for_spot(size(ordered_points, 1)*repeanting + 1 ) = total_time_for_spot(size(ordered_points, 1)*repeanting ) + btw_field_time;
                      cumulative_dose_at_point(size(ordered_points, 1)*repeanting + 1 ) = cumulative_dose_at_point(size(ordered_points, 1)*repeanting );
                  end
              end


              % Find the treshold
              upper_threshold = 1 - threshold/2;
              lower_threshold = threshold/2;

              % find the time within the treshold of the total dose
              target_dose_down_percent = lower_threshold * cumulative_dose_at_point(end);
              target_dose_up_percent = upper_threshold * cumulative_dose_at_point(end);
              cumulative_dose_at_point = round(cumulative_dose_at_point, 4);
              cumulative_dose_at_point = [0.0; cumulative_dose_at_point];
              total_time_for_spot = [ 0.0; total_time_for_spot];


              % Interpolate to get the exact time and dose
              [unique_dose, unique_indices] = unique(cumulative_dose_at_point);
              unique_total_time = total_time_for_spot(unique_indices);
              time_at_down_percent = interp1(unique_dose,  unique_total_time, target_dose_down_percent);
              time_at_up_percent = interp1(unique_dose, unique_total_time, target_dose_up_percent);
              dose_at_down_percent = interp1(unique_total_time, unique_dose, time_at_down_percent);
              dose_at_up_percent = interp1(unique_total_time, unique_dose, time_at_up_percent);

              time_interval = time_at_up_percent - time_at_down_percent;
              dose_interval = dose_at_up_percent - dose_at_down_percent;

              pbs_dr(point) = dose_interval/(0.001*time_interval);
              disp(['spot no :', num2str(point), ' -> pbs DR: ', num2str(pbs_dr(point)), ' Gy/s']);

              % Plot the cumulative dose at the point of interest as a function of time
              plot(total_time_for_spot, cumulative_dose_at_point, 'o-');

              hold on

          end

          disp(' ');
          pbs_dr_mean = mean(pbs_dr(ordered_points_in_target));
          disp(['mean pbs dose rate in target: ', num2str(mean(pbs_dr(ordered_points_in_target))), ' Gy/s']);
          disp(['0.5*(max+min) pbs dose rate: ', num2str(0.5*(max(pbs_dr)+min(pbs_dr))), ' Gy/s']);
          lgd = legend(string(2:1:size(point_of_interest,1)));
          legend('Location','northeastoutside','NumColumns',3)
          title(lgd,'spot no')
          title('Cumulative Dose at a Point as a Function of Time');
          xlabel('Time (ms)');
          ylabel('Cumulative Dose at Point of Interest (Gy)');
          grid on
          hold off

      end

      %% 3c) calcualte PBS dose rate at all point of interest

      if pbs_dose_rate == true

          x_range = linspace(0,(no_spots_X -1)*distance_X,no_spots_X);
          y_range = linspace(0,(no_spots_Y-1)*distance_Y,no_spots_Y);
          [XX, YY] = meshgrid(x_range(1)-1:0.1:x_range(end)+1, y_range(1)-1:0.1:y_range(end)+1);

          num_voxels = numel(XX);

          % Initialize dose rate matrix
          pbs_dr_grid = zeros(size(XX));

          % Iterate over all voxel positions
          for voxel_idx = 1:num_voxels
              [row, col] = ind2sub(size(XX), voxel_idx);
              voxel_position = [XX(row, col), YY(row, col)];

              cumulative_dose = zeros(size(ordered_points, 1) * no_repainting, 1);
              total_time_for_spot = zeros(size(ordered_points, 1) * no_repainting, 1);

              for repeanting = 1:no_repainting

                  if repeanting == 1
                      total_time_for_spot(1) = beam_on_spot * corner_weight * weight / no_repainting;
                      cumulative_dose(1) = corner_weight * weight / no_repainting * ...
                          exp(-((voxel_position(1) - ordered_points(1,1))^2 / (2 * sigma_X^2) + ...
                          (voxel_position(2) - ordered_points(1,2))^2 / (2 * sigma_Y^2)));
                  end

                  for i = 2:size(ordered_points, 1)
                      position = ordered_points(i-1, :);

                      % Determine weight factor
                      if (position(1) == 0 || position(1) == (no_spots_X-1)  * distance_X) && (position(2) == 0 || position(2) == (no_spots_Y-1)  * distance_Y)
                          weight_factor = corner_weight;
                      elseif (position(1) == 0 || position(1) == (no_spots_X-1)  * distance_X || position(2) == 0 || position(2) == (no_spots_Y-1)  * distance_Y)
                          weight_factor = rim_weight;
                      else
                          weight_factor = 1;
                      end

                      % Compute dose contribution at the voxel
                      dose_at_voxel = weight_factor * weight/no_repainting * ...
                          exp(-((voxel_position(1) - position(1))^2 / (2 * sigma_X^2) + ...
                          (voxel_position(2) - position(2))^2 / (2 * sigma_Y^2)));

                      % Accumulate dose values
                      idx = i + size(ordered_points, 1) * (repeanting - 1);
                      cumulative_dose(idx) = cumulative_dose(idx - 1) + dose_at_voxel;

                      % Calculate the time
                      switch_time_T = switch_X * (mod(i-1, (no_spots_X) ) > 0);
                      switch_time_U = switch_Y * (mod(i-1, (no_spots_Y) ) == 0 && i > 1);
                      total_time_for_spot(idx) = total_time_for_spot(idx - 1) + switch_time_T + switch_time_U + weight_factor * weight/no_repainting * beam_on_spot;
                  end

                  % Add time between repaintings
                  if no_repainting > 1 && no_repainting > repeanting
                      idx = size(ordered_points, 1) * repeanting + 1;
                      total_time_for_spot(idx) = total_time_for_spot(idx - 1) + btw_field_time;
                      cumulative_dose(idx) = cumulative_dose(idx - 1);
                  end
              end

              % Compute dose rate using thresholding
              upper_threshold = 1 - threshold/2;
              lower_threshold = threshold/2;

              target_dose_down_percent = lower_threshold * cumulative_dose(end);
              target_dose_up_percent = upper_threshold * cumulative_dose(end);
              cumulative_dose = round([0.0; cumulative_dose], 4);
              total_time_for_spot = [0.0; total_time_for_spot];

              % Interpolate time and dose values
              [unique_dose, unique_indices] = unique(cumulative_dose);
              unique_total_time = total_time_for_spot(unique_indices);
              time_at_down_percent = interp1(unique_dose, unique_total_time, target_dose_down_percent);
              time_at_up_percent = interp1(unique_dose, unique_total_time, target_dose_up_percent);
              dose_at_down_percent = interp1(unique_total_time, unique_dose, time_at_down_percent);
              dose_at_up_percent = interp1(unique_total_time, unique_dose, time_at_up_percent);

              % Compute dose rate for the voxel
              time_interval = time_at_up_percent - time_at_down_percent;
              dose_interval = dose_at_up_percent - dose_at_down_percent;

              pbs_dr_grid(row, col) = dose_interval / (0.001 * time_interval);
          end

          % Compute Mean Dose Rate in Isodose
          dose_on_dr_grid = interp2(xx, yy, dose_distribution, XX, YY, 'linear', 0);
          mask_level = dose_on_dr_grid >= level * max(dose_on_dr_grid(:));
          dose_rate_level = pbs_dr_grid(mask_level);
          mean_folk_dose_rate_90 = mean(dose_rate_level, 'omitnan');
          disp(['Mean Dose Rate in 90% Isodose Region: ', num2str(mean_folk_dose_rate_90), ' Gy/s']);

          % Plot the dose rate heatmap
          figure;
          mask = (XX >= 0 & XX <= (no_spots_X-1) * distance_X) & (YY >= 0 & YY <= (no_spots_Y-1) * distance_Y);
          pbs_dr_grid(~mask) = NaN; % Set outside region to NaN
          contourf(XX, YY, pbs_dr_grid, 200, 'LineColor', 'none');
          colormap('jet');
          hold on

          % Plot spot positions
          plot(ordered_points(:,1),ordered_points(:,2), ':k', 'MarkerEdgeColor', 'k');
          scatter(ordered_points(:,1),ordered_points(:,2), 'k', 'MarkerEdgeColor', 'k');
          scatter(ordered_points(1,1),ordered_points(1,2), 'rd', 'filled', 'MarkerEdgeColor', 'k');
          scatter(ordered_points(end,1),ordered_points(end,2), 'rx', 'filled', 'MarkerEdgeColor', 'k');
          scatter(ordered_points(18,1),ordered_points(18,2), '*', 'filled', 'MarkerEdgeColor', 'b');
          contour(xx,yy,dose_distribution, [level*max(max(dose_distribution)),level*max(max(dose_distribution))], 'LineColor', 'k');

          hold off;
          colorbar;
          grid on
          xlabel('X Position (cm)');
          ylabel('Y Position (cm)');
          title('ADR Distribution (Gy/s)');
          set(gca, 'YDir', 'normal');
      end

      %% 4) calcualte VdW dose rate at point of interest

      if vdw_dose_rate_spot_in_target == true

          cumulative_dose2R_at_point = zeros(size(ordered_points, 1)*no_repainting, 1);
          cumulative_dose_at_point = zeros(size(ordered_points, 1)*no_repainting, 1);
          point_of_interest = ordered_points;

          disp(' ');
          disp('---- VDW Average dose rate -----');

          for point = 1:size(ordered_points_in_target,2)
              for repeanting = 1:no_repainting

                  for i = 2:size(ordered_points, 1)

                      % calculate cumulative dose
                      position = ordered_points(i-1, :);
                      if (position(1) == 0 || position(1) == (no_spots_X-1)  * distance_X) && (position(2) == 0 || position(2) == (no_spots_Y-1)  * distance_Y)
                          weight_factor = corner_weight;
                      elseif position(1) == 0 || position(1) == (no_spots_X-1)  * distance_X || position(2) == 0 || position(2) == (no_spots_Y-1)  * distance_Y
                          weight_factor = rim_weight;
                      else
                          weight_factor = 1;
                      end

                      dose_at_point_of_interest = weight_factor * weight/no_repainting * exp(-((point_of_interest(point,1) - ordered_points(i-1, 1)).^2 / (2 * sigma_X^2) + (point_of_interest(point,2) - ordered_points(i-1, 2)).^2 / (2 * sigma_Y^2)));
                      cumulative_dose2R_at_point(i+ size(ordered_points, 1)*(repeanting-1)) = cumulative_dose2R_at_point(i+ size(ordered_points, 1)*(repeanting-1)-1) + dose_at_point_of_interest*dose_at_point_of_interest/( weight_factor * weight/no_repainting )*1/(beam_on_spot);
                      cumulative_dose_at_point(i+ size(ordered_points, 1)*(repeanting-1)) = cumulative_dose_at_point(i+ size(ordered_points, 1)*(repeanting-1)-1) + dose_at_point_of_interest;
                  end
              end

              dr(point) = cumulative_dose2R_at_point(end)/cumulative_dose_at_point(end)*1e3;

          end
          mean_dr = mean(dr);
          disp(['mean dose rate in target: ', num2str(mean(dr)), ' Gy/s']);
          disp(['0.5*(max+min) dose rate: ', num2str(0.5*(max(dr)+min(dr))), ' Gy/s']);

      end

      %% 4) calcualte VdW dose rate at all point of interest

      if vdw_dose_rate == true

          % Define voxel grid
          x_range = linspace(0,(no_spots_X -1)*distance_X,no_spots_X);
          y_range = linspace(0,(no_spots_Y-1)*distance_Y,no_spots_Y);
          [X, Y] = meshgrid(x_range(1)-1.5:0.05:x_range(end)+1.5, y_range(1)-1.5:0.05:y_range(end)+1.5);
          num_voxels = numel(X);

          % Initialize dose rate matrices
          cumulative_dose2R_grid = zeros(size(X));
          cumulative_dose_grid = zeros(size(X));
          dose_rate_grid = zeros(size(X));

          disp(' ');
          disp('---- VDW Average dose rate for entire grid -----');

          % Iterate over all voxel positions in the grid
          for voxel_idx = 1:num_voxels
              [row, col] = ind2sub(size(X), voxel_idx);
              voxel_position = [X(row, col), Y(row, col)];

              cumulative_dose2R = 0;
              cumulative_dose = 0;

              for repeanting = 1:no_repainting
                  for i = 2:size(ordered_points, 1)

                      % Get position of the dose spot
                      position = ordered_points(i-1, :);

                      % Determine weight factor
                      if (position(1) == 0 || position(1) == (no_spots_X-1)  * distance_X) && (position(2) == 0 || position(2) == (no_spots_Y-1)  * distance_Y)
                          weight_factor = corner_weight;
                      elseif position(1) == 0 || position(1) == (no_spots_X-1)  * distance_X || position(2) == 0 || position(2) == (no_spots_Y-1)  * distance_Y
                          weight_factor = rim_weight;
                      else
                          weight_factor = 1;
                      end

                      % Compute dose contribution at the voxel
                      dose_at_voxel = weight_factor * weight/no_repainting * exp(-((voxel_position(1) - position(1))^2 / (2 * sigma_X^2) + (voxel_position(2) - position(2))^2 / (2 * sigma_Y^2)));

                      % Accumulate dose values
                      cumulative_dose2R = cumulative_dose2R + dose_at_voxel^2 / ( weight_factor *weight/no_repainting) * 1/(beam_on_spot);
                      cumulative_dose = cumulative_dose + dose_at_voxel;
                  end
              end

              % Compute dose rate for the voxel
              dose_rate_grid(row, col) = (cumulative_dose2R / cumulative_dose) * 1e3;
          end

          % Compute and display mean dose rate
          mean_dr = mean(dose_rate_grid(:), 'omitnan');
          disp(['Mean dose rate in target: ', num2str(mean_dr), ' Gy/s']);
          disp(['0.5*(max+min) dose rate: ', num2str(0.5 * (max(dose_rate_grid(:)) + min(dose_rate_grid(:)))), ' Gy/s']);

          % Compute Mean Dose Rate in Isodose
          dose_on_dr_grid = interp2(xx, yy, dose_distribution, X, Y, 'linear', 0);
          mask_level = dose_on_dr_grid >= level * max(dose_on_dr_grid(:));
          dose_rate_level = dose_rate_grid(mask_level);
          mean_vdw_dose_rate_level = mean(dose_rate_level, 'omitnan');
          disp(['Mean DVW Dose Rate in 90% Isodose Region: ', num2str(mean_vdw_dose_rate_level), ' Gy/s']);

          % Plot the dose rate heatmap
          figure;
          mask = (X >= 0 & X <= (no_spots_X-1) * distance_X) & (Y >= 0 & Y <= (no_spots_Y-1) * distance_Y);
          dose_rate_grid(~mask) = NaN; % Set outside region to NaN
          contourf(X, Y, dose_rate_grid, 100, 'LineColor', 'none');
          colormap('jet');
          hold on

          % Plot spot positions
          plot(ordered_points(:,1),ordered_points(:,2), ':k', 'MarkerEdgeColor', 'k');
          scatter(ordered_points(:,1),ordered_points(:,2), 'k', 'MarkerEdgeColor', 'k');
          scatter(ordered_points(1,1),ordered_points(1,2), 'rd', 'filled', 'MarkerEdgeColor', 'k');
          scatter(ordered_points(end,1),ordered_points(end,2), 'rx', 'filled', 'MarkerEdgeColor', 'k');

          % Plot point of interest
          scatter(ordered_points(18,1),ordered_points(18,2), '*', 'filled', 'MarkerEdgeColor', 'b');

          % Plot rectangle at 95% dose position
          contour(xx,yy,dose_distribution, [level*max(max(dose_distribution)),level*max(max(dose_distribution))], 'LineColor', 'k');

          hold off;
          colorbar;
          grid on
          xlabel('X Position (cm)');
          ylabel('Y Position (cm)');
          title('DADR Distribution (Gy/s)');
          set(gca, 'YDir', 'normal');

      end

      %% 5) calculate the field dose rate

      if field_dose_rate == true

          dose_to_consider =  max(max(dose_distribution)); %can be changed
          time = ( switch_X*(no_spots_X-1)  + switch_Y*(no_spots_Y-1)  + dose_to_consider/no_repainting * beam_on_spot * ((no_spots_Y-1) +(no_spots_X-1) ) )* no_repainting + btw_field_time * (no_repainting - 1);
          dead_time = ( switch_X*(no_spots_X-1)  + switch_Y*(no_spots_Y-1)  )* no_repainting + btw_field_time * (no_repainting - 1);
          beam_on = time-dead_time;
          adr = dose_to_consider/(0.001*time);

          disp(' ');
          disp('---- Average dose rate -----');
          disp(['total irradiation time: ', num2str(time), ' ms']);
          disp(['total dead time time: ', num2str(dead_time), ' ms']);
          disp(['total beam on time: ', num2str(beam_on), ' ms']);
          disp(['total dose: ', num2str(dose_to_consider), ' Gy']);
          disp(['average dose rate: ', num2str(adr), ' Gy/s']);

      end
end
