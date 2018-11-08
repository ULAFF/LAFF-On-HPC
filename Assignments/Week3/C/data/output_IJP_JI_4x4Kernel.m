Name = 'IJP\_JI\_4x4Kernel';
 
% number of repeats:% 3
% enter first, last, inc:% 48 960 48 
data = [
%  n          reference      |         current implementation 
%        time       GFLOPS   |    time       GFLOPS     diff 
   960 3.0613e-01 5.7801e+00    1.0194e+00 1.7358e+00 9.9476e-13
   912 1.5271e-01 9.9344e+00    2.8616e-01 5.3016e+00 8.8107e-13
   864 7.4199e-02 1.7385e+01    1.5667e-01 8.2338e+00 8.5265e-13
   816 9.4843e-02 1.1458e+01    2.5712e-01 4.2264e+00 7.6739e-13
   768 5.5649e-02 1.6280e+01    1.5914e-01 5.6929e+00 7.6739e-13
   720 5.8473e-02 1.2766e+01    1.0647e-01 7.0111e+00 7.1054e-13
   672 2.9384e-02 2.0655e+01    1.0203e-01 5.9486e+00 6.5370e-13
   624 2.2495e-02 2.1602e+01    5.5756e-02 8.7156e+00 5.6843e-13
   576 1.9367e-02 1.9735e+01    7.0415e-02 5.4279e+00 4.8317e-13
   528 3.8331e-02 7.6803e+00    1.2040e-01 2.4451e+00 3.9790e-13
   480 6.4457e-02 3.4315e+00    1.1785e-01 1.8768e+00 3.6948e-13
   432 1.7975e-02 8.9704e+00    1.0854e-01 1.4855e+00 3.1264e-13
   384 8.6659e-03 1.3068e+01    8.9785e-02 1.2613e+00 2.2737e-13
   336 8.7222e-03 8.6980e+00    2.1916e-02 3.4617e+00 1.7053e-13
   288 7.8291e-03 6.1024e+00    8.6384e-03 5.5306e+00 1.1369e-13
   240 1.8198e-03 1.5193e+01    5.4239e-03 5.0974e+00 4.2633e-14
   192 8.1468e-04 1.7376e+01    1.7813e-03 7.9471e+00 2.8422e-14
   144 4.0394e-04 1.4784e+01    1.7692e-03 3.3754e+00 2.8422e-14
    96 1.4304e-04 1.2370e+01    2.0204e-04 8.7581e+00 1.0658e-14
    48 3.5822e-05 6.1745e+00    2.6518e-05 8.3409e+00 7.1054e-15
];

% Maximum difference between reference and your implementation: 9.947598e-13.