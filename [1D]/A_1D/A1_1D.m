classdef A1_1D
    methods (Static)
        %% > 1. -----------------------------------------------------------
        % >> 1.1. ---------------------------------------------------------
        %  > Set up input variables #1.
        function [inp] = Set_msh(el,t)
            inp.h   = el;     %  > Edge length.
            inp.Lim = [-1,1]; %  > Mesh limits (x-direction).
            inp.t   = t;      %  > Mesh type.
            if inp.t == 2 || inp.t == 3
                switch inp.t
                    case 2, inp.x = [5,0];
                    case 3, inp.x = [1];
                end
            end
        end
        % >> 1.2. ---------------------------------------------------------
        %  > Set up input variables #2.
        function [inp] = Set_inp(x)
            %  > ----------------------------------------------------------
            %  > Boundary conditions.
            inp.B(1)  = "Dirichlet";   %  > West(W).
            inp.B(2)  = "Dirichlet";   %  > East(E).
            if ~all(ismember(inp.B,["Dirichlet","Neumann","Robin"]))
                return;
            end
            %  > ----------------------------------------------------------
            %  > Coefficients/function handle(s).
            inp.C     = A1_1D.Fh_1;
            inp.F.Fh  = A1_1D.Fh_2(x);
            %  > ----------------------------------------------------------
            %  > Test #.
            inp.T     = 1;
            %  > ----------------------------------------------------------
            %  > Problem setup.
            inp.PS.EE = 1;
            inp.PS.P  = 1;
            %  > ----------------------------------------------------------
        end
        %  > 1.2.1. -------------------------------------------------------
        function [F] = Fh_1()
            C = [1;-1];
            F = [cell(size(C))];
            for i = 1:size(C,1)
                for j = 1:size(C,2)
                    F{i,j} = @(x) repmat(C(i),size(x));
                end
            end
        end
        %  > 1.2.2. -------------------------------------------------------
        function [F] = Fh_2(x)
            switch x
                case 1
                    a = 5;
                    F = @(x) x(:,1).^a;
                case 2
                    a = [1e+02];
                    b = [1e+02];
                    c = [0];
                    f = [cell(size(c))];
                    for i = 1:numel(c)
                        f{i} = @(x) a(i).*exp(-b(i).*(x(:,1)-c(i)).^2);
                    end
                    F = @(x) sum(cellfun(@(f) f(x),f));
                case 3
                    a = [1e+02];
                    b = [1e+02];
                    c = [0];
                    d = [-pi./8,+pi./8];
                    f = [cell(size(c))];
                    g = [cell(size(d))];
                    for i = 1:numel(c)
                        f{i} = @(x) a(i).*exp(-b(i).*(x(:,1)-c(i)).^2);
                    end
                    for i = 1:numel(d)
                        g{i} = @(x) sqrt((x(:,1)-d(i)).^2);
                    end
                    F = @(x) sum(cellfun(@(f) f(x),cat(2,f,g)));
                otherwise
                    return;
            end
        end
    end
end