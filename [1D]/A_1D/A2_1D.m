classdef A2_1D
    methods (Static)
        %% > 1. -----------------------------------------------------------
        % >> 1.1. ---------------------------------------------------------
         %  > Set up 'msh' structure.
        function [msh] = Set_msh(el,t)
            %  > ----------------------------------------------------------
            % >> Xv.
            Xv = A2_1D.Set_Xv(A1_1D.Set_msh(el,t));
            %  > ----------------------------------------------------------
            % >> c.
            msh.c.Nc = numel(Xv)-1;
            for i = 1:msh.c.Nc
                msh.c.f.if  (i,:) = [i,i+1];
                msh.c.f.Nf  (i,:) = [-1,1];
                msh.c.Xc    (i,1) = 1./2.*(Xv(i)+Xv(i+1));
                msh.c.Volume(i,1) = Xv(i+1)-Xv(i);
            end
            %  > ----------------------------------------------------------
            % >> f.
            msh.f.Nf = numel(Xv);
            msh.f.l  = true (msh.f.Nf,1);
            for i = 1:msh.f.Nf
                msh.f.ic{i,1} = [i-1,i];
                if i == 1 || i == msh.f.Nf
                    msh.f.ic{i,1} = setdiff(msh.f.ic{i,1},[0,msh.f.Nf]);
                    msh.f.l (i,1) = false;
                end
            end
            msh.f.Xv     = Xv';
            msh.f.Volume = ones(msh.f.Nf,1);
            %  > ----------------------------------------------------------
            % >> h.
            msh.h = 1./msh.c.Nc;
            %  > ----------------------------------------------------------
            %  > Sort fields...
            msh = A2_1D.sort(msh);
        end
        
        %% > 2. -----------------------------------------------------------
        % >> 2.1. ---------------------------------------------------------
        %  > switch inp.t
        %  >    case 1, uniformly spaced mesh.
        %  >    case 2, smooth non-uniform mesh.
        %  >    case 3, highly distorted non-uniform mesh.
        %  > end
        function [Xv] = Set_Xv(inp)
            %  > Auxiliary variables.
            dX = diff (inp.Lim);
            Nc = round(dX./inp.h);
            ni = 0:1/Nc:1;
            Xi = inp.Lim(1)+dX.*ni;
            
            %  > Select mesh type...
            switch inp.t
                case 1
                    Xv      = Xi;
                case 2
                    dx      = inp.x(2)-inp.Lim(1);
                    i (1)   = 1-(1-exp( inp.x(1))).*dx./dX;
                    i (2)   = 1-(1-exp(-inp.x(1))).*dx./dX;
                    B       = 1./(2.*inp.x(1)).*log(i(1)./i(2));
                    j (1,:) = inp.x(1).*(ni-B);
                    j (2,:) = inp.x(1).*B;
                    Xv      = inp.Lim(1)+dx.*(1+sinh(j(1,:))./sinh(j(2,:)));
                case 3
                    rng default;
                    i       = [1,Nc+1];
                    j       = AuxiliaryFunctions.setdiff(1:Nc+1,i);
                    Xv(i)   = Xi(i);
                    Xv(j)   = Xi(j)+(rand(1,Nc-1)-0.5).*(1./Nc).^inp.x(1);
            end
        end

        %% > 3. -----------------------------------------------------------
        % >> 3.1. ---------------------------------------------------------
        %  > Sort 'msh' structure fields.
        function [msh] = sort(msh)
            msh     = orderfields(msh    ,{'c','f','h'});
            msh.c   = orderfields(msh.c  ,{'f','Nc','Volume','Xc'});
            msh.c.f = orderfields(msh.c.f,{'if','Nf'});
            msh.f   = orderfields(msh.f  ,{'ic','l','Nf','Volume','Xv'});
        end
    end
end