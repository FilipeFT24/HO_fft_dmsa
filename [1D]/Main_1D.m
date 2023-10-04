%% > Clear memory, clean screen, close any figure.
clear, clc, close all; beep off; warning on;
%% > Run...
% >> ----------------------------------------------------------------------
%  > Working directories.
addpath(genpath('A_1D'));
addpath(genpath('B_1D'));
addpath(genpath( '../'));
% >> ----------------------------------------------------------------------
el  = 1.00e-01; %  > Edge length.
mt  = 1;        %  > Mesh type.
ms  = 1;        %  > Manufactured solution (MMS).
inp = A1_1D.Set_inp(ms);
obj = B2_1D.Run_P  (inp,el,mt);
% >> ----------------------------------------------------------------------