%% Performance of implementations of loop orderings with all outer loops
%% This Live Script
% We now look at the different implementations around calls to the BLAS routines 
% "dgemv" (matrix-vector multiplication) and "dger" (rank-1 update).
% 
% To gather the performance data, in the command (terminal) window change 
% the directory to LAFF-On-HPC/Assignments/Week1/C/.  After implementing the various 
% versions,  execute
% 
%       make Plot_All_Outer
% 
% (This reruns many previous experiments, but with larger problem sizes.)

plot_colors = [ 0 0 0; 0 0 1; 0 1 0; 0 1 1; 1 0 0; 1 0 1; 1 1 0; 1 1 1];

% Create figure
figure1 = figure('Name','GFLOPS');

% Create axes, labels, legends.  In future routines for plotting performance, 
% the next few lines will be hidden in the script.
axes2 = axes('Parent',figure1);
hold(axes2,'on');
ylabel( 'GFLOPS', 'FontName', 'Helvetica Neue' );
xlabel( 'matrix dimension m=n=k', 'FontName', 'Helvetica Neue' );
box(axes2,'on');
set( axes2, 'FontName', 'Helvetica Neue', 'FontSize', 18);
             
% Plot time data for I_bli_dgemv 
output_I_bli_dgemv   % load data for JPI ordering
assert( max(abs(data(:,6))) < 1.0e-10, ...
    'Hmmm, better check if there is an accuracy problem');
plot( data(:,1), data(:,5), 'DisplayName', 'I\_bli\_dgemv ', 'MarkerSize', 8, 'LineWidth', 2, ...
      'Marker', 'o', 'LineStyle', '-.', 'Color', plot_colors( 2,: ) );

% Plot time data for I_dgemv (to plot change "0" to "1")
if ( 0 ) 
  output_I_dgemv
  assert( max(abs(data(:,6))) < 1.0e-10, ...
      'Hmmm, better check if there is an accuracy problem');
  plot( data(:,1), data(:,5), 'DisplayName', 'I\_dgemv', 'MarkerSize', 8, 'LineWidth', 2, ...
        'Marker', 'o', 'LineStyle', '-.', 'Color', plot_colors( 3,: ) );
end

% Plot time data for J_bli_dgemv  (to plot change "0" to "1")
if ( 1 ) 
  output_J_bli_dgemv
  assert( max(abs(data(:,6))) < 1.0e-10, ...
      'Hmmm, better check if there is an accuracy problem');
  plot( data(:,1), data(:,5), 'DisplayName', 'J\_bli\_dgemv', 'MarkerSize', 8, 'LineWidth', 2, ...
        'Marker', 'o', 'LineStyle', '-', 'Color', plot_colors( 4,: ) );
end

% Plot time data for J_dgemv  (to plot change "0" to "1")
if ( 0 ) 
  output_J_dgemv
  assert( max(abs(data(:,6))) < 1.0e-10, ...
      'Hmmm, better check if there is an accuracy problem');
  plot( data(:,1), data(:,5), 'DisplayName', 'J\_dgemv', 'MarkerSize', 8, 'LineWidth', 2, ...
        'Marker', 'o', 'LineStyle', '-', 'Color', plot_colors( 5,: ) );
end

% Plot time data for P_bli_dger (to plot change "0" to "1")
if ( 1 ) 
  output_P_bli_dger
  assert( max(abs(data(:,6))) < 1.0e-10, ...
      'Hmmm, better check if there is an accuracy problem');
  plot( data(:,1), data(:,5), 'DisplayName', 'P\_bli\_dger', 'MarkerSize', 8, 'LineWidth', 3, ...
        'Marker', 'o', 'LineStyle', '--', 'Color', plot_colors( 6,: ) );
end

% Plot time data for P_dger  (to plot change "0" to "1")
if ( 0 ) 
  output_P_dger
  assert( max(abs(data(:,6))) < 1.0e-10, ...
      'Hmmm, better check if there is an accuracy problem');
  plot( data(:,1), data(:,5), 'DisplayName', 'P\_dger', 'MarkerSize', 8, 'LineWidth', 3, ...
        'Marker', 'o', 'LineStyle', '--', 'Color', plot_colors( 7,: ) );
end

% Optionally show the reference implementation performance data
if ( 0 )
  plot( data(:,1), data(:,3), 'MarkerSize', 8, 'LineWidth', 1, ...
        'LineStyle', '-', 'DisplayName', 'Ref', 'Color', plot_colors( 1,: ) );
end

% Adjust the x-axis and y-axis range to start at 0
v = axis;                   % extract the current ranges
axis( [ 0 v(2) 0 v(4) ] )   % start the x axis and y axis at zero

% Optionally change the top of the graph to capture the theoretical peak
if ( 0 )
    turbo_clock_rate = 2.6;
    flops_per_cycle = 16;
    peak_gflops = turbo_clock_rate * flops_per_cycle;

    axis( [ 0 v(2) 0 peak_gflops ] )  
end

legend2 = legend( axes2, 'show' );
set( legend2, 'Location', 'northeast', 'FontSize', 18) ;

% Uncomment if you want to create a pdf for the graph
% print( 'Plot_Inner_P.pdf', '-dpdf' );
%%