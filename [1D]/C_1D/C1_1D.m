%% > Clear memory, clean screen, close any figure.
clear, clc, close all; beep off; warning on;
%% > Run...
% >> ----------------------------------------------------------------------
%  > Working directories.
addpath(genpath( '../'));
% >> ----------------------------------------------------------------------
ms = 1; %  > Manufactured solution (MMS).
T  = 1; %  > Test #.
switch T
    case 1
        %  > 'inp'.
        inp = A1_1D.Set_inp(1);
        %  > 'msh' and 'obj'.
        a        = 2;
        b        = 0.5;
        msh      = A2_1D.Set_msh(2.00E-02,1);
        n        = round(b.*msh.f.Nf);
        obj(a+1) = struct();
        for i = 1:a+1
            switch i
                case 1
                    obj(i).c = 0;
                    obj(i).d = 0;
                    obj(i).f = A3_1D.Initialize_f(inp,msh);
                    obj(i).s = B1_1D.Initialize_s(inp,msh);
                    obj(i).s = B1_1D.Update_s    (inp,msh,obj(i).f,obj(i).s);
                    obj(i).e = NE_XD.Initialize_e(inp,msh);
                    obj(i).e = NE_XD.Update_e    (    msh,obj(i).e,obj(i).s);
                otherwise
                    switch i
                        case 2
                            x = sort(obj(1).e.a.t.f_abs(:,1),'descend'); c = x(n+1);
                            d = find(obj(1).e.a.t.f_abs(:,1) >= c);
                        case 3
                            x = sort(obj(1).e.a.t.f_abs(:,3),'descend'); c = x(n);
                            d = find(obj(1).e.a.t.f_abs(:,3) >= c);
                        otherwise
                            return;
                    end
                    obj(i)   = obj(1);
                    obj(i).c = c;
                    obj(i).d = d;
                    for j = reshape(setdiff(string(fieldnames(obj(i).s)),"u"),1,[])
                        for k = 1:numel(obj(i).s.u.p.(j))
                            obj(i).s.u.p.(j){k}(d,:) = 2.*ceil(obj(i).s.u.p.(j){k}(d,:)./2)+1;
                            obj(i).s.u.s.(j){k}      = d;
                        end
                    end
                    obj(i).s = B1_1D.Update_s(inp,msh,obj(i).f,obj(i).s);
                    obj(i).e = NE_XD.Update_e(    msh,obj(i).e,obj(i).s);
            end
        end
        %  > Plot...
        A  = cell(size(obj,2)-1);
        for i = 1:size(obj,2)
            if i ~= 1
                A{i-1} = msh.f.Xv([obj(i).d([true;diff(obj(i).d) ~= 1]),obj(i).d([diff(obj(i).d) ~= 1;true])]);
            end
        end
        Z1 = {{A{1},A{2}},{A{1},{}},{{},A{2}}};
        for i = 1:size(obj,2)
            PP_1D.Plot(inp,msh,obj(i),...
                {[0,0,0,0,0],...
                 [1,1],...
                 [0],...
                 [0,0],...
                 [0,0],...
                 [1],...
                 [0],...
                 [0,0,0,0]},Z1{i});
        end
        Z2 = cell(size(obj,2)-1,1);
        for i = 1:size(obj,2)-1
            Z2{i} = [obj(i+1).e.a.n_abs.t.f(:,1)./obj(1).e.a.n_abs.t.f(:,1),...
                     obj(i+1).e.a.n_abs.t.f(:,3)./obj(1).e.a.n_abs.t.f(:,3),...
                     obj(i+1).e.a.n_abs.t.c(:,1)./obj(1).e.a.n_abs.t.c(:,1),...
                     obj(i+1).e.a.n_abs.  c(:,1)./obj(1).e.a.n_abs.  c(:,1)];
        end
        PP_1D.Plot(inp,msh,obj(end),...
            {[0,0,0,0,0],...
             [0,0],...
             [1],...
             [1,1],...
             [0,0],...
             [0],...
             [0],...
             [0,0,0,0]},Z2);
    case 2
        %  > 'inp'.
        inp       = A1_1D.Set_inp(ms);
        inp.PS.EE = 1;
        %  > 'msh' and 'obj'.
        el        = exp(1).^(linspace(log(1.00e-02),log(1.00e-03),5));
        for i = 1:numel(el)
            msh(i) = A2_1D.Set_msh(el(i),1); %#ok<*SAGROW>
            for j  = 1:numel(inp)
                obj(i,j).f = A3_1D.Initialize_f(inp(j),msh(i));
                obj(i,j).s = B1_1D.Initialize_s(inp(j),msh(i));
                obj(i,j).s = B1_1D.Update_s    (inp(j),msh(i),obj(i,j).f,obj(i,j).s);
                obj(i,j).e = NE_XD.Initialize_e(inp(j),msh(i));
                obj(i,j).e = NE_XD.Update_e    (       msh(i),obj(i,j).e,obj(i,j).s);
                V  (i,j).e = obj(i,j).e;
                V  (i,j).h = msh(i)  .h;
            end
        end
        %  > Plot...
        PP_1D.Plot(inp,msh,V,...
            {[0,0,0,0,0],...
             [0,0],...
             [0],...
             [0,0],...
             [0,0],...
             [0],...
             [0],...
             [1,1,1,1]},{[]});
    otherwise
        return;
end
% >> ----------------------------------------------------------------------