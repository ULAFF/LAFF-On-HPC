%% Performance of various blocked matrix-matrix multiplications.
% We next look at the different implementations that block the operands into 
% submatrices, with a "micro-kernel"  that is progressively further optimized.
% 
% To gather the performance data, in the command (terminal) window change 
% the directory to Assignments/Week2/C/.  After implementing a new implementation, 
% execute it with "make" to gather data.  You will want to activate the plotting 
% of a new performance curve by changing "( 0 )" to "( 1 )" in the appropriate 
% place.

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

system( 'cp ../../../Week1/C/data/output_IJP.m .' );
output_IJP   % load data for IJP ordering
assert( max(abs(data(:,6))) < 1.0e-10, ...
    'Hmmm, better check if there is an accuracy problem');
plot( data(:,1), data(:,5), 'DisplayName', 'IJP', 'MarkerSize', 8, 'LineWidth', 2, ...
    'Marker', 'o', 'LineStyle', '-.', 'Color', plot_colors( 2,: ) );
  
% Plot performance data for IJP_mbxnbxkb
if ( 1 )
  output_IJP_mbxnbxkb   % load data for JPI ordering
  assert( max(abs(data(:,6))) < 1.0e-10, ...
      'Hmmm, better check if there is an accuracy problem');
  plot( data(:,1), data(:,5), 'DisplayName', 'IJP\_mbxnbxkb', 'MarkerSize', 8, 'LineWidth', 2, ...
      'Marker', 'o', 'LineStyle', '-.', 'Color', plot_colors( 3,: ) );
end

% Plot time data for IJ_mbxnbxk (to plot change "0" to "1")
if ( 0 ) 
  output_IJ_mbxnbxk
  assert( max(abs(data(:,6))) < 1.0e-10, ...
      'Hmmm, better check if there is an accuracy problem');
  plot( data(:,1), data(:,5), 'DisplayName', 'IJ\_mbxnbxk', 'MarkerSize', 8, 'LineWidth', 2, ...
        'Marker', 'o', 'LineStyle', '-.', 'Color', plot_colors( 4,: ) );
end

% Plot time data for JI_mbxnbxk  (to plot change "0" to "1")
if ( 0 ) 
  output_JI_mbxnbxk
  assert( max(abs(data(:,6))) < 1.0e-10, ...
      'Hmmm, better check if there is an accuracy problem');
  plot( data(:,1), data(:,5), 'DisplayName', 'JI\_mbxnbxk', 'MarkerSize', 8, 'LineWidth', 2, ...
        'Marker', 'o', 'LineStyle', '-', 'Color', plot_colors( 5,: ) );
end

% Plot time data for IJ_P_bli_dger (to plot change "0" to "1")
if ( 1 ) 
  output_IJ_P_bli_dger
  assert( max(abs(data(:,6))) < 1.0e-10, ...
      'Hmmm, better check if there is an accuracy problem');
  plot( data(:,1), data(:,5), 'DisplayName', 'IJ\_P\_bli\_dger', 'MarkerSize', 8, 'LineWidth', 3, ...
        'Marker', 'o', 'LineStyle', '--', 'Color', plot_colors( 6,: ) );
end

% Plot time data for IJ_4x4Kernel  (to plot change "0" to "1")
if ( 0 ) 
  output_IJ_4x4Kernel
  assert( max(abs(data(:,6))) < 1.0e-10, ...
      'Hmmm, better check if there is an accuracy problem');
  plot( data(:,1), data(:,5), 'DisplayName', 'IJ\_4x4Kernel', 'MarkerSize', 8, 'LineWidth', 3, ...
        'Marker', 'o', 'LineStyle', '--', 'Color', plot_colors( 7,: ) );
end

% Plot time data for JI_4x4Kernel  (to plot change "0" to "1")
if ( 0 ) 
  output_JI_4x4Kernel
  assert( max(abs(data(:,6))) < 1.0e-10, ...
      'Hmmm, better check if there is an accuracy problem');
  plot( data(:,1), data(:,5), 'DisplayName', 'JI\_4x4Kernel', 'MarkerSize', 8, 'LineWidth', 3, ...
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
% print( 'Plot_Blocked_MMM.pdf', '-dpdf' );
%%