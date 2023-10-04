classdef A3_1D
    methods (Static)
        %% > 1. -----------------------------------------------------------
        % >> 1.1. ---------------------------------------------------------
        %  > Initialize field 'f'.
        function [f] = Initialize_f(inp,msh)
            %  > 'fh' (function handles).
            f.fh = A3_1D.Update_fh(inp);
            %  > 'bd' (boundary face indices/surface normals/type).
            f.bd = A3_1D.Update_bd(inp,msh,f.fh);
            %  > 'st' ((volumetric) source term).
            f.st = A3_1D.Update_st(msh,f.fh.func);
        end
        %  > 1.1.1. -------------------------------------------------------
        %  > Update field 'fh' (function handles).
        function [fh] = Update_fh(inp)
            %  > Symbolic variables.
            x       = sym('x');
            f       = inp.F.Fh(x);
            %  > Convert to function handle...
            df      = diff(f,x);
            fh.d    = matlabFunction(df,'Vars',x);
            fh.f    = matlabFunction(f ,'Vars',x);
            %  > 'func'.
            fh.func = inp.C{1}(x).*f+inp.C{2}(x).*df;
        end
        %  > 1.1.2. -------------------------------------------------------
        %  > Update field 'bd' (boundary face indices/surface normals/type).
        function [bd] = Update_bd(inp,msh,fh)
            %  > 'i'.
            bd.i = cat(1,1,msh.f.Nf);
            %  > 'Nf'.
            for i = 1:size(bd.i,1)
                bd.Nf(i,1) = msh.c.f.Nf(msh.f.ic{bd.i(i,1),1},msh.c.f.if(msh.f.ic{bd.i(i,1),1},:) == bd.i(i,1),:);
            end
            %  > 't' and 'v'.
            bd.t = inp.B';
            for i = 1:size(bd.t,1)
                switch bd.t(i,1)
                    case "Dirichlet"
                        bd.v(i,1) = fh.f(msh.f.Xv(bd.i(i,1)));
                    case {"Neumann","Robin"}
                        switch bd.t(i,1)
                            case "Neumann"
                                V3 = fh.d(msh.f.Xv(bd.i(i,1)));
                            case "Robin"
                                V1 = fh.f(msh.f.Xv(bd.i(i,1)));
                                V2 = inp.C{2}(msh.f.Xv(bd.i(i,1),1))./inp.C{1}(msh.f.Xv(bd.i(i,1),1)).*fh.d(msh.f.Xv(bd.i(i,1),1));
                                V3 = V1+V2;
                        end
                        bd.v(i,1) = V3*bd.Nf(i,1);
                end
            end
        end
        %  > 1.1.3. -------------------------------------------------------
        %  > Update field 'st' ((volumetric) source term).
        function [st] = Update_st(msh,func)
            %  > Symbolic variables.
            x  = sym('x');
            %  > Convert to function handle...
            Fh = matlabFunction(func,'Vars',x);
            %  > For each cell...
            st = zeros(msh.c.Nc,size(msh.c.f.if,2));
            for i = 1:msh.c.Nc
                for j = 1:size(msh.c.f.if,2)
                    st(i,j) = Fh(msh.f.Xv(msh.c.f.if(i,j),1)).*msh.c.f.Nf(i,j);
                end
            end
            st(isinf(st) | isnan(st)) = 0;
        end
    end
end