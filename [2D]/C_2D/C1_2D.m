%% > Clear memory, clean screen, close any figure.
clear, clc, close all; beep off; warning on;
%% > Run...
% >> ----------------------------------------------------------------------
%  > Working directories.
addpath(genpath( '../'));
% >> ----------------------------------------------------------------------
wf = 1; %  > Weighting function.
ms = 2; %  > Manufactured solution (MMS).
T  = 1; %  > Test #.
switch T
    case 1
        %  > Auxiliary variables.
        LD(1) = 1;
        LD(2) = 1;
        MS    = "[2D]/C_2D/[.mat Files]/[.T1]/[.M]/M1_2D.mat"; MD = dir(MS);
        VS    = "[2D]/C_2D/[.mat Files]/[.T1]/[.V]/A1_2D.mat"; VD = dir(VS);
        %  > 'msh'.
        if ~LD(1)
            el             =    exp(1).^(linspace(log(4.00e-02),log(2.00e-02),5));
            msh(numel(el)) = struct();
            for i = 1:numel(el)
                msh(i) = A2_2D.Set_msh(el(i),1); fprintf("Cycle #%3d (msh)\n",i);
            end
            save(MS,"msh");
            M2 = msh;
        else
            M1 = struct2cell(load(MD.name));
            M2 = cat(2,M1{:});
        end
        %  > 'inp'.
        inp       = A1_2D.Set_inp(ms,wf);
        inp.PS.EE = 1;
        %  > 'obj'.
        if ~LD(2)
            obj(size(M2,2),numel(inp)) = struct();
            V  (size(M2,2),numel(inp)) = struct();
            for i = 1:size(M2,2)
                for j = 1:numel(inp)
                    obj(i,j).f = A3_2D.Initialize_f(inp(j),M2(i));
                    obj(i,j).s = B1_2D.Initialize_s(inp(j),M2(i));
                    obj(i,j).s = B1_2D.Update_s    (inp(j),M2(i),obj(i,j).f,obj(i,j).s);
                    obj(i,j).e = NE_XD.Initialize_e(inp(j),M2(i));
                    obj(i,j).e = NE_XD.Update_e    (       M2(i),obj(i,j).e,obj(i,j).s);
                    V  (i,j).e = obj(i,j).e;
                    V  (i,j).h = M2 (i)  .h;
                    fprintf("Cycle #%3d (obj)\n",i);
                end
            end
            save(VD,"V");
            V2 = V;
        else
            V1 = struct2cell(load(VD.name));
            V2 = cat(2,V1{:});
        end
        %  > Plot...
        PP_2D.Plot(inp,M2,V2,...
            {[0,0],...
             [0,0],...
             [0,0],...
             [0],...
             [0,0],...
             [1,1,1]});
    otherwise
        return;
end
% >> ----------------------------------------------------------------------