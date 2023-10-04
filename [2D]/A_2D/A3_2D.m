classdef A3_2D
    methods (Static)
        %% > 1. -----------------------------------------------------------
        % >> 1.1. ---------------------------------------------------------
        %  > Initialize field 'f'.
        function [f] = Initialize_f(inp,msh)
            %  > 'fh' (function handles).
            f.fh = A3_2D.Update_fh(inp);
            %  > 'bd' (boundary face indices/surface normals/type).
            f.bd = A3_2D.Update_bd(inp,msh,f.fh);
            %  > 'qd' (1D quadrature).
            f.qd = @(u,x) x(1,:).*(1-u)./2+x(2,:).*(1+u)./2;
            %  > 'st' ((volumetric) source term).
            f.st = A3_2D.Update_st(msh,f.fh.func);
        end
        %  > 1.1.1. -------------------------------------------------------
        %  > Update field 'fh' (function handles).
        function [fh] = Update_fh(inp)
            %  > Symbolic variables.
            x          = sym('x');
            y          = sym('y');
            f          = inp.F.Fh([x,y]);
            %  > Convert to function handle...
            df     {1} = diff(f,x);
            df     {2} = diff(f,y);
            fh.d   {1} = matlabFunction(df{1},'Vars',{[x,y]});
            fh.d   {2} = matlabFunction(df{2},'Vars',{[x,y]});
            fh.f       = matlabFunction(f    ,'Vars',{[x,y]});
            fh.func{1} = inp.C{1,1}([x,y]).*f+inp.C{2,1}([x,y]).*df{1};
            fh.func{2} = inp.C{1,2}([x,y]).*f+inp.C{2,2}([x,y]).*df{2};
        end
        %  > 1.1.2. -------------------------------------------------------
        %  > Update field 'bd' (boundary face indices/surface normals/type).
        function [bd] = Update_bd(inp,msh,fh)
            %  > 'i'.
            bd.i = find(~msh.f.l);
            %  > 'Nf'.
            for i = 1:sum(~msh.f.l,1)
                bd.Nf(i,:) = msh.c.f.Nf{msh.f.ic{bd.i(i,1),1},1}(msh.c.f.if(msh.f.ic{bd.i(i,1),1},:) == bd.i(i,1),:);
            end
            %  > 't' and 'v'.
            for i = 1:sum(~msh.f.l,1)
                %  > 't'.
                n = bd.Nf(i,:);
                if     n(1,2) > -n(1,1) && n(1,2) <  n(1,1), bd.t(i,1) = inp.B(1); %  > East (E).
                elseif n(1,2) >  n(1,1) && n(1,2) > -n(1,1), bd.t(i,1) = inp.B(2); %  > North(N).
                elseif n(1,2) < -n(1,1) && n(1,2) >  n(1,1), bd.t(i,1) = inp.B(3); %  > West (W).
                elseif n(1,2) <  n(1,1) && n(1,2) < -n(1,1), bd.t(i,1) = inp.B(4); %  > South(S).
                else
                    return;
                end
                %  > 'v'.
                xc = msh.f.xy.c(bd.i(i,1),:);
                switch bd.t(i,1)
                    case "Dirichlet"
                        bd.v(i,1) = fh.f(xc);
                    case {"Neumann","Robin"}
                        switch bd.t(i,1)
                            case "Neumann"
                                V3 = zeros(1,size(n,2));
                                for j = 1:size(n,2)
                                    V3(1,j) = fh.d{j}(xc);
                                end
                            case "Robin"
                                V1 = repmat(fh.f(xc),1,size(n,2));
                                V2 = zeros (1,size(n,2));
                                for j = 1:size(inp.C,2)
                                    V2(1,j) = inp.C{2,j}(xc)./inp.C{1,j}(xc).*fh.d{j}(xc);
                                end
                                V2(isinf(V2) | isnan(V2)) = 0; V3 = V1+V2;
                        end
                        bd.v(i,1) = V3*bd.Nf(i,:)';
                end
            end
        end
        %  > 1.1.3. -------------------------------------------------------
        %  > Update field 'st' ((volumetric) source term).
        function [st] = Update_st(msh,F)
            %  > Symbolic variables.
            c   = sym('c',[2,size(msh.f.xy.c,2)]);
            t   = sym('t');
            x   = sym('x');
            y   = sym('y');
            %  > r(t).
            r   = (1-t).*(c(1,:))+t.*(c(2,:));
            %  > F(r(t)).
            Ft  = [subs(-F{2},{x,y},{r(1),r(2)}),...
                   subs( F{1},{x,y},{r(1),r(2)})];
            %  > F(r(t))*dr(t).
            Fdr = dot(Ft,diff(r,t,1));
            %  > Convert to function handle...
            Fh  = matlabFunction(Fdr,'Vars',{c,t});
            Fh  = @(c)(@(t) Fh(c,t));
            %  > For each cell...
            st  = zeros(msh.c.Nc,size(msh.struct.ConnectivityList,2));
            for i = 1:msh.c.Nc
                for j = 1:size(msh.c.c.xy.v{i},1)
                    if j ~= size(msh.c.c.xy.v{i},1)
                        k = [j,j+1];
                    else
                        k = [j,1];
                    end
                    st(i,j) = integral(Fh(msh.c.c.xy.v{i}(k,:)),0,1);
                end
            end
            st(isinf(st) | isnan(st)) = 0;
        end
        %  > 1.1.4. -------------------------------------------------------
        %  > 1D quadrature.
        function [Q] = Q_1D(p)
            %  > From 'quadGaussLegendre(n,varargin)'...
            n     = max (ceil(p./2),[],2);
            B     = sqrt(((1:n-1)./(2:n))./((2*(0:n-2)+1)./((0:n-2)+1).*(2*(1:n-1)+1)./((1:n-1)+1)));
            J     = diag(B,1)+diag(B,-1);
            %  > Abcissas/weights.
            [V,D] = eig (J,'vector');
            [X,I] = sort(D,1); 
            %  > Assign to structure 'Q'.
            Q.x   = X';
            Q.w   = V(1,I).^2;
        end
    end
end