classdef NE_XD
    methods (Static)
        %% > 1. -----------------------------------------------------------
        % >> 1.1. ---------------------------------------------------------
        %  > Initialize field 'e'.
        function [e] = Initialize_e(inp,msh)
            %  > Auxiliary variables.
            a = msh.c.Nc; %  > Nc.
            b = msh.f.Nf; %  > Nf.
            c = 2;        %  > Number of norms (L1 and L\infty).
            if ~inp.PS.EE
                fn = ["a"];
            else
                fn = ["a","d","x"];
            end
            
            %  > Loop through...
            for i = fn
                switch i
                    case "a", x    = 1; % >> w/ analytic solution.
                              y    = 3; %  > a0,a1*h^n,a0+a1*h^n+...
                    case "d", x    = 3; % >> difference between...
                              y(1) = 3; %  > a0,a1*h^n,a0+a1*h^n+... %  > Φ   -φ(l).
                              y(2) = 3; %  > a0,a1*h^n,a0+a1*h^n+... %  > Φ   -φ(h).
                              y(3) = 3; %  > a0,a1*h^n,a0+a1*h^n+... %  > φ(h)-φ(l).
                    case "x", x    = 2; % >> w/ discrete solutions.
                              y(1) = 1; %  > a1*h^n.
                              y(2) = 3; %  > a0,a1*h^n,a0+a1*h^n+...
                end
                for j = 1:x
                    %  > Error distribution.
                    e.(i)(j).c.c       = zeros(a,1);
                    e.(i)(j).c.c_abs   = zeros(a,1);
                    e.(i)(j).t.c       = zeros(a,1);
                    e.(i)(j).t.c_abs   = zeros(a,1);
                    e.(i)(j).t.f_abs   = zeros(b,y(j));
                    %  > Error norms.
                    e.(i)(j).n_abs.c   = zeros(c,1);
                    e.(i)(j).n_abs.t.c = zeros(c,1);
                    e.(i)(j).n_abs.t.f = zeros(c,y(j));
                end
            end
        end
        % >> 1.2. ---------------------------------------------------------
        %  > Update field 'e'.
        %  > 1.2.1. -------------------------------------------------------
        %  > Update error distribution/norms.
        function [e] = Update_e(msh,e,s)
            %  > Loop through...
            for i = reshape(string(fieldnames(e)),1,[])
                switch i
                    case "a"
                        dx.(i){1} = s.c.x.x.dx.a;
                    case "d"
                        dx.(i){1} = s.c.x.x.dx.a   -s.c.x.x.dx.x{1}; %  > Φ   -φ(l).
                        dx.(i){2} = s.c.x.x.dx.a   -s.c.x.x.dx.x{2}; %  > Φ   -φ(h).
                        dx.(i){3} = s.c.x.x.dx.x{2}-s.c.x.x.dx.x{1}; %  > φ(h)-φ(l).
                    case "x"
                        dx.(i)    = s.c.x.x.dx.x;
                end
                for j = 1:numel(dx.(i))
                    %  > Truncation error distribution.
                    e.(i)(j).t       = NE_XD.ed(msh,dx.(i){j},e.(i)(j).t);
                    %  > Discretization error distribution.
                    e.(i)(j).c.c     = s.c.m.A\e.(i)(j).t.c;
                    e.(i)(j).c.c_abs = abs(e.(i)(j).c.c);
                    %  > Error norms.
                    e.(i)(j)         = NE_XD.en(msh,e.(i)(j));
                end
            end
        end
        %  > 1.2.2. -------------------------------------------------------
        %  > Auxiliary function #1.
        function [e] = ed(msh,dx,e)
            %  > c.
            e.c     = sum(dx(:,:,1),2);
            e.c_abs = abs(e.c);
            %  > f.
            for i = 1:msh.f.Nf
                j            = 1;
                k            = msh.f.ic{i,1}(j);
                e.f_abs(i,:) = abs(dx(k,msh.c.f.if(k,:) == i,:));
            end
        end
        %  > 1.2.3. -------------------------------------------------------
        %  > Auxiliary function #2.
        function [e] = en(msh,e)
            %  > c.
            e.n_abs.c   = [NE_XD.Lk(e.c.c_abs,1,msh.c.Volume);max(e.c.c_abs,[],1)];
            e.n_abs.t.c = [NE_XD.Lk(e.t.c_abs,1,msh.c.Volume);max(e.t.c_abs,[],1)];
            %  > f.
            e.n_abs.t.f = [NE_XD.Lk(e.t.f_abs,1,msh.f.Volume);max(e.t.f_abs,[],1)];
        end
        %  > 1.2.4. -------------------------------------------------------
        %  > Auxiliary function #3.
        function [L] = Lk(e,k,A)
            L = (sum(A.*e.^k,1)./sum(A,1)).^(1./k);
        end
    end
end