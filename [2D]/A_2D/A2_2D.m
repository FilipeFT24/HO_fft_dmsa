classdef A2_2D
    methods (Static)
        %% > 1. -----------------------------------------------------------
        % >> 1.1. ---------------------------------------------------------
        %  > Set up 'msh' structure.
        function [msh] = Set_msh(el)
            % >> struct.
            msh.struct     = A2_2D.struct_2D(el);
            % >> c.
            CL_c           = msh.struct.ConnectivityList;
            msh.c.Nc       = size(CL_c,1);
            %  > c.c.
            %  #1: Cell 'i': face/vertex neighbours + auxiliary array 'V'.
            [msh.c.c.nb]   = A2_2D.cc_nb   (CL_c);
            %  #2: Cell 'i': centroid/vertex coordinates (for each cell).
            msh.c.c.xy     = A2_2D.cc_xy   (msh.struct);
            %  #3: Cell 'i': volume.
            msh.c.Volume   = A2_2D.c_Volume(msh.c.c.xy.v);
            %  > c.f.
            %  #1: Identify all/boundary/bulk faces.
            F{1}           = A2_2D.cf_F1   (CL_c,msh.c.c.nb.f);
            %  #2: Cell 'i': outer face normals.
            msh.c.f        = A2_2D.cf_NSf  (F{1},msh.struct,msh.c.c.xy.c);
            % >> f.
            %  #1: List faces.
            F{2}           = A2_2D.f_F2(F{1});
            %  -1) Number of (unique) faces.
            msh.f.Nf       = size(F{2}.ic,1);
            %  -2) Face 'i': cell indices.
            msh.f.ic       = F{2}.ic;
            %  -3) Face 'i': vertex indices.
            msh.f.iv       = F{2}.iv;
            %  -4) Face 'i': is it a boundary or bulk face?
            msh.f.l        = F{2}.l;
            %  #2: Face 'i': centroid/vertex coordinates.
            msh.f.xy       = A2_2D.f_xy (msh.f.iv,msh.struct);
            %  #3: Cell 'i': face indices.
            msh.c.f.if     = A2_2D.cf_if(F,msh.struct);
            %  #4: Face 'i': volume.
            msh.f.Volume   = A2_2D.f_Volume(msh.f.xy.v);
            %  #5: Face 'i': min. distance.
            aux = [msh.c.c.xy.c;msh.f.xy.c(~msh.f.l,:)];

            aa = round(sqrt((msh.f.xy.c(:,1)-aux(:,1)').^2+(msh.f.xy.c(:,2)-aux(:,2)').^2),round(-log(el),0)+1);
            bb(1,1:msh.c.Nc) = 1:msh.c.Nc;
            bb(1,1+msh.c.Nc:msh.c.Nc+sum(~msh.f.l,1)) = find(~msh.f.l);
            bb(2,1:msh.c.Nc) = 1;
            bb(2,1+msh.c.Nc:msh.c.Nc+sum(~msh.f.l,1)) = 0;

            msh.d = {aa,bb};
            
            % >> d.
            %  #1: Domain edge/reference length.
            msh.h          = 1./sqrt(msh.c.Nc);
            %  > Sort fields...
            msh            = A2_2D.sort(msh);
        end
        
        %% > 2. -----------------------------------------------------------
        % >> 2.1. ---------------------------------------------------------
        function [struct] = struct_2D(el)
            %  > inp.
            inp = A1_2D.Set_msh(el);
            %  > Select cell polyhedral (type) and construct mesh...
            if inp.p == 1
                %  > XYv.
                [V.X,V.Y] = meshgrid(linspace(inp.Lim(1,1),inp.Lim(1,2),round(diff(inp.Lim(1,:),1,2)./inp.el)+1),...
                                     linspace(inp.Lim(2,1),inp.Lim(2,2),round(diff(inp.Lim(2,:),1,2)./inp.el)+1));
                %  > Connectivity list and points.
                Nc(1) = size(V.X,2)-1;
                Nc(2) = size(V.Y,1)-1;
                for i = 1:Nc(1)
                    for j = 1:Nc(2)
                        struct.ConnectivityList(Nc(2).*(i-1)+j,1) = (Nc(2)+1).*(i-1)+j;   % > SW.
                        struct.ConnectivityList(Nc(2).*(i-1)+j,2) = (Nc(2)+1).*(i-0)+j;   % > SE.
                        struct.ConnectivityList(Nc(2).*(i-1)+j,3) = (Nc(2)+1).*(i-0)+j+1; % > NE.
                        struct.ConnectivityList(Nc(2).*(i-1)+j,4) = (Nc(2)+1).*(i-1)+j+1; % > NW.
                    end
                end
                struct.Points = cat(2,reshape(V.X,[],1),reshape(V.Y,[],1));
            else
                %  > Auxiliary variables.
                dLim = diff (inp.Lim,1,2);
                ncX  = round(dLim(1)./el,0);
                ncY  = round(2*ncX/sqrt(3),0);
                nvX  = ncX+1;
                nvY  = ncY+1;

                %  > Initialize...
                struct.ConnectivityList = zeros(2*ncX*ncY,3);
                struct.Points           = zeros(round(nvY/2,0)*nvX+floor(nvY/2)*(nvX+1),2);
                %  > Connectivity list.
                for i = 1:ncY
                    switch rem(i-1,2)
                        case 0
                            A1 = (i-1)/2*(2*nvX+1);
                            A2 = (i-1)/2*(2*nvX+1);
                            B1 = (1:ncX);
                            B2 = (1:nvX);

                        case 1
                            A1 = (i-2)/2*(2*nvX+1)+nvX;
                            A2 = (i-2)/2*(2*nvX+1)+nvX+1;
                            B1 = (1:nvX);
                            B2 = (1:ncX);
                    end
                    struct.ConnectivityList((i-1)*(2*nvX-1)+1:i*(2*nvX-1),:) = [A1+[B1;B1+1;B1+nvX+1]';A2+[B2;B2+nvX+1;B2+nvX]'];
                end
                %  > Points.
                for i = 1:nvY
                    r = rem (i-1,2);
                    j = nvX*(i-1)+round(i/2,0)-1+(1:nvX+r);
                    switch r
                        case 0
                            struct.Points(j,1) = inp.Lim(1)+dLim(1)*(0:ncX)/ncX;
                        case 1
                            struct.Points(j,1) = inp.Lim(1)+dLim(1)*[0,0.5:ncX-0.5,ncX]/ncX;
                    end
                    struct.Points(j,2) = inp.Lim(2)+dLim(2)*repmat(i-1,nvX+r,1)/ncY;
                end
                %[struct.Points,struct.ConnectivityList] = ...
                %    distmesh2d(@(p) drectangle(p,inp.Lim(1,1),inp.Lim(1,2),inp.Lim(2,1),inp.Lim(2,2)),@(p) ones(size(p,1),1),inp.el,inp.Lim',[inp.Lim(:,1)';diag(inp.Lim')';diag(fliplr(inp.Lim),0)';inp.Lim(:,2)'],[1E-03,1E-03]);
            end
        end
        
        %% > 3. -----------------------------------------------------------
        % >> 3.1. ---------------------------------------------------------
        %  > Field: 'c'.
        %  > 3.1.1. -------------------------------------------------------
        %  > #1: Cell 'i': face/vertex neighbours+auxiliary array 'V'.
        function [nb] = cc_nb(CL_c)
            %  > Auxiliary variables.
            CL_s    = sort (CL_c,2);
            sz      = size (CL_c);
            sz  (3) = max  (CL_c,[],"all");
            f_ij    = cell (sz(1),1);
            V   {1} = false(sz(1),sz(3));
            V   {2} = false(sz(1),sz(1));
            
            %  > Assemble sparse matrix...
            a     = repelem(1:sz(1),sz(2))';
            b     = reshape(CL_c',[],1);
            V{1}  = sparse (a,b,true(size(a,1),1));
            %  > Vertex neighbours.
            l     = 1:sz(1);
            for i = 1:sz(1)
                for j = 1:sz(2)
                    V{2}(i,V{1}(:,CL_c(i,j))) = true;
                end
                V     {2}(i,i) = false;
                nb.v{i,1}(:,1) = l(V{2}(i,:));
            end
            %  > Face neighbours.
            for i = 1:sz(1)
                j = size(nb.v{i,1},1);
                for k = 1:j
                    f_ij{i,1}(k,1) = sum(ismembc(CL_s(i,:),CL_s(nb.v{i}(k),:)),2) > 1;
                end
                nb.f{i,1}(:,1) = nb.v{i,1}(f_ij{i,1});
            end
        end
        %  > 3.1.2. -------------------------------------------------------
        %  > #2: Cell 'i': centroid/vertex coordinates (for each cell).
        function [xy] = cc_xy(struct)
            for i = 1:size(struct.ConnectivityList,1)
                j              = 1:size(struct.Points,2);
                xy.v{i,1}(:,j) = struct.Points(struct.ConnectivityList(i,:),j);
                xy.c     (i,j) = AuxiliaryFunctions.mean (xy.v{i,1}(:,j),1);
            end
        end
        %  > 3.1.3. -------------------------------------------------------
        %  > #3: Cell 'i': volume.
        function [Volume] = c_Volume(xy_v)
            Volume = zeros(size(xy_v,1),1);
            for i = 1:size(xy_v,1)
                Volume(i,1) = polyarea(xy_v{i}(:,1),xy_v{i}(:,2));
            end
        end
        %  > 3.1.4. -------------------------------------------------------
        %  > #3: Cell 'i': volume.
        function [Volume] = f_Volume(xy_v)
            Volume = zeros(size(xy_v,1),1);
            for i = 1:size(xy_v,1)
                Volume(i,1) = sqrt(diff(xy_v{i}(:,1),[],1).^2+diff(xy_v{i}(:,2),[],1).^2);
            end
        end
        %  > %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %  > 3.1.5. -------------------------------------------------------
        %  > Identify all/boundary/bulk faces (only face neighbours are evaluated).
        function [F] = cf_F1(CL_c,nb_f)
            %  > Auxiliary variables.
            CL_s = sort(CL_c,2);
            sz   = size(CL_c);
            fb   = cell(sz(1),1);
            fi   = cell(sz(1),1);
            fv   = cell(sz(1),1);
            vj   = cat (1,1:sz(2),circshift(1:sz(2),sz(2)-1))';
            n    = size(vj,2);

            for i = 1:sz(1)
                %  > Deal all faces (through connectivity matrix).
                for j = 1:sz(2)
                    fv{i,1}(j,1:n) = sort(CL_c(i,vj(j,:)),2);
                    fv{i,1}(j,n+1) = i;
                end
                fi{i,1} = false(sz(2),1);
                %  > Check and increment accordingly (if bulk face)...
                l = 1;
                for j = 1:size(nb_f{i,1},1)
                    k = ismembc(CL_s(i,:),CL_s(nb_f{i,1}(j),:));
                    %  > If the evaluated cells have more than 1 vertex in common...
                    if sum(k,2) > 1
                        fb{i,1}(l) = find(all(bsxfun(@eq,fv{i,1}(:,1:end-1),CL_s(i,k)),2)); l = l+1;
                    end
                end
                %  > fi: 0-bnd.
                %        1-blk.
                fi{i,1}(fb{i,1},1) = true;
            end
            %  > Assign to structure 'F'.
            F.c.i = fi;
            F.c.v = fv;
            F.f.i = cat(1,F.c.i{:,1});
            F.f.v = cat(1,F.c.v{:,1});
        end
        %  > 3.1.6. -------------------------------------------------------
        %  > #2: Cell 'i': outer face normals (Sf).
        function [f] = cf_NSf(F1,struct,xy_cm)
            xy_v = cell (size(F1.c.v,1),1);
            v_Sf = zeros(size(struct.ConnectivityList,2),2);
            for i = 1:size(F1.c.v,1)
                %  > xy_v.
                for j = 1:size(F1.c.v{i,1},1)
                    xy_v{i,1}{j,1} = struct.Points(F1.c.v{i}(j,1:end-1),:);
                end
                %  > Sf.
                for j = 1:size(xy_v{i,1},1)
                    %  > \vec{FC}.
                    v_FC      = xy_cm(i,:)-AuxiliaryFunctions.mean(xy_v{i,1}{j,1},1);
                    %  > \vec{Sf}.
                    v_Sf(j,2) = xy_v{i,1}{j,1}(1,1)-xy_v{i,1}{j,1}(2,1);
                    v_Sf(j,1) = xy_v{i,1}{j,1}(2,2)-xy_v{i,1}{j,1}(1,2);
                    %  > Check...
                    if dot(v_FC,v_Sf(j,:)) > 0
                        v_Sf(j,:) = -v_Sf(j,:);
                    end
                end
                f.Nf{i,1} = bsxfun(@rdivide,v_Sf,sqrt(sum(v_Sf.^2,2)));
                f.Sf{i,1} = v_Sf;
            end
        end
        %  > %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % >> 3.2. ---------------------------------------------------------
        %  > Field: 'f'.
        %  > 3.2.1. -------------------------------------------------------
        %  > Auxiliary function (identify vertex indices of bulk faces).
        function [blk] = blk_f(F1)
            %  > Auxiliary variables.
            f_blk = F1.f.v(F1.f.i,:);
            
            %  > Face 'i' belongs to cells 'j(1)' and 'j(2)'.
            [~,~,a] = unique(f_blk(:,1:end-1),'rows');
            [~,b,~] = unique(a);
            for i = 1:size(b,1)
                blk.if(i,:) = find(a == a(b(i,1),1),2)';
            end
            blk.if = sortrows(blk.if,1);
            %  > Get indices...
            for i = 1:size(blk.if,1)
                blk.ic(i,:) = f_blk(blk.if(i,:),end);
            end
        end
        %  > 3.2.2. -------------------------------------------------------
        %  > List faces.
        function [F2] = f_F2(F1)
            %  > Auxiliary variables.
            fb = A2_2D.blk_f(F1);
            l  = 1;
            n  = size(fb.if,1);
            
            %  > List of unique faces.
            j     = find(F1.f.i);
            k     = 1;
            [a,b] = sortrows(cat(1,j(fb.if(:,k)),find(~F1.f.i)),1);
            %  > Assign to structure 'F2'.
            F2.iv           = F1.f.v(a(:,1),1:end-1);
            F2.l (b >  n,1) = false;
            F2.l (b <= n,1) = true;
            for i = 1:size(a,1)
                if ~F2.l (i,1)
                    F2.ic{i,1} = F1.f.v(a(i,1),end);
                else
                    F2.ic{i,1} = fb.ic(l,:); l = l+1;
                end
            end
        end
        %  > %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %  > 3.2.3. -------------------------------------------------------
        %  > #2: Face 'i': centroid/vertex coordinates.
        function [xy] = f_xy(f,struct)
            for i = 1:size(f,1)
                j              = 1:size(f,2);
                xy.v{i,1}(j,:) = struct.Points(f(i,j),:);
                xy.c     (i,j) = AuxiliaryFunctions.mean(xy.v{i,1}(:,j),1);
            end
        end
        %  > %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %  > 3.2.4. -------------------------------------------------------
        %  > #3: Cell 'i': face indices.
        function [ic_f] = cf_if(F,struct)
            %  > Auxiliary variables.
            sz_c = size (F{1}.c.i,1);
            sz_f = size (F{2}.ic ,1);
            ic_f = zeros(size(struct.ConnectivityList));
            U    = false(sz_c,sz_f);
            j    = 1:sz_f;
            
            for i = j
                U(F{2}.ic{i,1},i) = true;
            end
            for i = 1:sz_c
                ic_f(i,:) = j(U(i,:));
                c         = A2_2D.AB(F{2}.iv(ic_f(i,:),:),F{1}.c.v{i,1}(:,1:end-1));
                ic_f(i,c) = ic_f(i,:);
            end
        end
        %  > 3.2.4.1. -----------------------------------------------------
        %  > Auxiliary function (faster version of built-in function 'intersect' w/ 'rows' option).
        function [l] = AB(A,B)
            [~,i] = sortrows(A);
            [~,j] = sortrows(B);
            [~,k] = sortrows(i);
              [l] = j(k);
        end
        
        %% > 4. -----------------------------------------------------------
        % >> 4.1. ---------------------------------------------------------
        %  > Sort 'msh' structure fields.
        function [msh] = sort(msh)
            msh        = orderfields(msh       ,{'c','d','f','h','struct'});
            msh.c      = orderfields(msh.c     ,{'c','f','Nc','Volume'});
            msh.c.c    = orderfields(msh.c.c   ,{'nb','xy'});
            msh.c.c.nb = orderfields(msh.c.c.nb,{'f','v'});
            msh.c.c.xy = orderfields(msh.c.c.xy,{'c','v'});
            msh.c.f    = orderfields(msh.c.f   ,{'if','Nf','Sf'});
            msh.f      = orderfields(msh.f     ,{'ic','iv','l','Nf','Volume','xy'});
            msh.f.xy   = orderfields(msh.f.xy  ,{'c','v'});
        end
    end
end