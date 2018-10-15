%% Performance of implementations of loop orderings with outer loop indexed with "P"
%% This Live Script
% This Live Script helps you visualize the performance of implementations that 
% order the loops so that the "J" loop is the outer-most loop:  Gemm_PIJ.c, Gemm_PJI.c, 
% Gemm_P_Ger_J_Axpy, Gemm_P_Ger_I_Axpy, etc.
% 
% To gather the performance data, in the command (terminal) window change 
% the directory to Assignments/Week1/C/.  After implementing the various versions,  
% execute 
% 
%        make PIJ   (actually, you probably did this one already)
% 
%        make PJI
% 
%        make P_Ger_J_Axpy
% 
%         make P_Ger_I_Axpy
% 
%        make P_dger
% 
%        make P_bli_dger
% 
% or, alternatively,
% 
%       make Plot_Outer_P
% 
% This compiles and executes a driver routine (the source of which is in 
% driver.c) that collects accuracy and performance data for the various implementations.  
% 
% When completed, various data is in output file 'output_XYZ.m' (for XYZ 
% $$ \in $$ {PIJ, PJI, ...}) in the same directory where you found this Live Script 
% (LAFF-On-HPC/Assignments/Week1/C/data/).  This Life Script then creates graphs 
% from that timing data.  Go ahead and click on "Run All".  It executes all the 
% code in the rest of this file.  You will want to look at the graphs this creates.

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
             
% Plot time data for PIJ 
output_PIJ   % load data for JPI ordering
assert( max(abs(data(:,6))) < 1.0e-10, ...
    'Hmmm, better check if there is an accuracy problem');
plot( data(:,1), data(:,5), 'DisplayName', 'PIJ', 'MarkerSize', 8, 'LineWidth', 2, ...
      'Marker', 'o', 'LineStyle', '-.', 'Color', plot_colors( 2,: ) );

% Plot time data for PJI (to plot change "0" to "1")
if ( 0 ) 
  output_PJI 
  assert( max(abs(data(:,6))) < 1.0e-10, ...
      'Hmmm, better check if there is an accuracy problem');
  plot( data(:,1), data(:,5), 'DisplayName', 'PJI', 'MarkerSize', 8, 'LineWidth', 2, ...
        'Marker', 'o', 'LineStyle', '-.', 'Color', plot_colors( 3,: ) );
end

% Plot time data for P_Ger_J_Axpy  (to plot change "0" to "1")
if ( 0 ) 
  output_P_Ger_J_Axpy
  assert( max(abs(data(:,6))) < 1.0e-10, ...
      'Hmmm, better check if there is an accuracy problem');
  plot( data(:,1), data(:,5), 'DisplayName', 'P\_Ger\_J\_Axpy', 'MarkerSize', 8, 'LineWidth', 2, ...
        'Marker', 'o', 'LineStyle', '-.', 'Color', plot_colors( 4,: ) );
end

% Plot time data for P_Ger_I_Axpy  (to plot change "0" to "1")
if ( 0 ) 
  output_P_Ger_I_Axpy
  assert( max(abs(data(:,6))) < 1.0e-10, ...
      'Hmmm, better check if there is an accuracy problem');
  plot( data(:,1), data(:,5), 'DisplayName', 'P\_Ger\_I\_Axpy', 'MarkerSize', 8, 'LineWidth', 2, ...
        'Marker', 'o', 'LineStyle', '-.', 'Color', plot_colors( 5,: ) );
end

% Plot time data for P_dger (to plot change "0" to "1")
if ( 0 ) 
  output_P_dger
  assert( max(abs(data(:,6))) < 1.0e-10, ...
      'Hmmm, better check if there is an accuracy problem');
  plot( data(:,1), data(:,5), 'DisplayName', 'P\_dger', 'MarkerSize', 8, 'LineWidth', 3, ...
        'Marker', 'o', 'LineStyle', '--', 'Color', plot_colors( 6,: ) );
end

% Plot time data for P_bli_dger  (to plot change "0" to "1")
if ( 0 ) 
  output_P_bli_dger
  assert( max(abs(data(:,6))) < 1.0e-10, ...
      'Hmmm, better check if there is an accuracy problem');
  plot( data(:,1), data(:,5), 'DisplayName', 'P\_bli\_dger', 'MarkerSize', 8, 'LineWidth', 3, ...
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
set( legend2, 'Location', 'northwest', 'FontSize', 18) ;0
% Uncomment if you want to create a pdf for the graph
% print( 'Plot_Outer_P.pdf', '-dpdf' );
%%