classdef A1_2D
    methods (Static)
        %% > 1. -----------------------------------------------------------
        % >> 1.1. ---------------------------------------------------------
        %  > Set up input variables #1.
        function [inp] = Set_msh(el)
            inp.el  = el;          %  > Edge length.
            inp.Lim = [-1,1;-1,1]; %  > Mesh limits (x/y-direction).
            inp.p   = 1;           %  > Cell polyhedral.
            if inp.p ~= 1 && inp.p ~= 2
                return;
            end
        end
        % >> 1.2. ---------------------------------------------------------
        %  > Set up input variables #2.
        function [inp] = Set_inp(x,y)
            %  > ----------------------------------------------------------
            %  > Boundary conditions.
            inp.B(1)  = "Dirichlet"; %  > East (E).
            inp.B(2)  = "Dirichlet"; %  > North(N).
            inp.B(3)  = "Dirichlet"; %  > West (W).
            inp.B(4)  = "Dirichlet"; %  > South(S).
            if ~all(ismember(inp.B,["Dirichlet","Neumann","Robin"]))
                return;
            end
            %  > ----------------------------------------------------------
            %  > Coefficients/function handle(s).
            inp.C     = A1_2D.Fh_1;
            inp.F.Fh  = A1_2D.Fh_2(x);
            inp.F.Wf  = A1_2D.Fh_3(y);
            %  > ----------------------------------------------------------
            %  > Test #.
            inp.T     = 1;
            %  > ----------------------------------------------------------
            %  > Problem setup.
            inp.PS.EE = 0;
            inp.PS.P  = 1;
            %  > ----------------------------------------------------------
        end
        %  > 1.2.1. -------------------------------------------------------
        function [F] = Fh_1()
            C = [[1,1];  %  > u.
                -[1,1]]; %  > Gamma.
            F = [cell(size(C))];
            for i = 1:size(C,1)
                for j = 1:size(C,2)
                    F{i,j} = @(x) repmat(C(i,j),size(x,1),1);
                end
            end
        end
        %  > 1.2.2. -------------------------------------------------------
        function [F] = Fh_2(x)
            switch x
                case 1
                    a = repmat( 5,1,2);
                    % = repmat(10,1,2);
                    F = @(x) x(:,1).^a(1).*x(:,2).^a(2);
                case 2
                    a = [1e+00];
                    b = [1e+02];
                    c = [0,0];
                    f = [cell(size(c,1))];
                    for i = 1:size(c,1)
                        f{i} = @(x) a(i).*exp(-b(i).*((x(:,1)-c(i,1)).^2+(x(:,2)-c(i,2)).^2));
                    end
                    F = @(x) sum(cellfun(@(f) f(x),f));
                case 3
                    a = [1e+00];
                    b = [1e+02];
                    c = [0,0];
                    d = [0,+pi./8;+pi./8,0];
                    f = [cell(size(c,1),1)];
                    g = [cell(size(d,1),1)];
                    for i = 1:size(c,1)
                        f{i} = @(x) a(i).*exp(-b(i).*((x(:,1)-c(i,1)).^2+(x(:,2)-c(i,2)).^2));
                    end
                    for i = 1:size(d,1)
                        g{i} = @(x) sqrt((x(:,1)-d(i,1)).^2+(x(:,2)-d(i,2)).^2);
                    end
                    F = @(x) sum(cellfun(@(f) f(x),cat(1,f,g)));
                otherwise
                    return;
            end
        end
        %  > 1.2.3. -------------------------------------------------------
        function [F] = Fh_3(x)
            switch x
                case 1
                    F = @(dx,n) sqrt(sum(dx.^2,2)).^(-n);
                otherwise
                    return;
            end
        end
    end
end