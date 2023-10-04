%% > Clear memory, clean screen, close any figure.
clear, clc, close all; beep off; warning on;
%% > Run...
% >> ----------------------------------------------------------------------
%  > Working directories.
addpath(genpath('A_2D'));
addpath(genpath('B_2D'));
addpath(genpath( '../'));
% >> ----------------------------------------------------------------------
el  = 1.00e-01; %  > Edge length.
ms  = 2;        %  > Manufactured solution (MMS).
wf  = 1;        %  > Weighting function.
inp = A1_2D.Set_inp( ms,wf);
obj = B2_2D.Run_P  (inp,el);
% >> ----------------------------------------------------------------------