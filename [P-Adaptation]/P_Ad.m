classdef P_Ad
    methods (Static)
        %% > 1. -----------------------------------------------------------
        % >> 1.1. ---------------------------------------------------------
        %  > 'P-Adaptive' run.
        function [OBJ2] = P_Adaptive(inp,msh)
            %  > Auxiliary variables/initialize...
            u = false(1,3);
            if size(inp.C,2) == 1
                LD = 0;
            else
                LD = 1;
                SD = "[P-Adaptation]/[.mat Files]/[2D]/[.T1]/A1_2D.mat";
            end
            inp.PS.EE     = 1;
            inp.PS.P      = 1;
            inp.P_Ad.ec   = 1.00e-10; %  > Min. discretization error (L_inf norm).
            inp.P_Ad.fNNZ = 2.85;     %  > Max. NNZ/NNZ_0.
            inp.P_Ad.fp   = 0.10;     %  > Fixed-fraction parameter.
            inp.P_Ad.nc   = 20;       %  > Max. number of cycles.
            inp.P_Ad.R2   = 1;        %  > Apply rule #2?
            u(1)          = 1;        %  > Φ.
            u(2)          = 1;        %  > φ(h).
            u(3)          = 1;        %  > φ(l).
            v             = find(u);
            
            if ~LD
                %  > 'obj'.
                switch size(inp.C,2)
                    case 1
                        obj.f = A3_1D.Initialize_f(inp,msh);
                        obj.s = B1_1D.Initialize_s(inp,msh);
                        obj.s = B1_1D.Update_s    (inp,msh,obj.f,obj.s);
                    case 2
                        obj.f = A3_2D.Initialize_f(inp,msh);
                        obj.s = B1_2D.Initialize_s(inp,msh);
                        obj.s = B1_2D.Update_s    (inp,msh,obj.f,obj.s);
                end
                obj.e = NE_XD.Initialize_e(inp,msh);
                obj.e = NE_XD.Update_e    (    msh,obj.e,obj.s);
                %  > Until any stopping criterion has been met...
                objI  = cell(1,sum(u,2));
                OBJI  = cell(1,sum(u,2));
                for i = 1:sum(u,2)
                    j                       = 1;
                    objI{i}                 = obj;
                    OBJI{i}(j).e            = objI{i}(j).e;
                    OBJI{i}(j).s.c.m.nnz.At = objI{i}(j).s.c.m.nnz.At;
                    OBJI{i}(j).s.u.p.c      = objI{i}(j).s.u.p.c;
                    while 1
                        %  > Check stopping criteria and update/save cycle count/solution...
                        if ~P_Ad.Check_SC(inp,objI{i},v(i))
                            j                       = j+1;
                            objI{i}(j)              = P_Ad.Update_P(inp,msh,objI{i}(j-1).e,objI{i}(j-1).f,objI{i}(j-1).s,v(i));
                            OBJI{i}(j).e            = objI{i}(j).e;
                            OBJI{i}(j).s.c.m.nnz.At = objI{i}(j).s.c.m.nnz.At;
                            OBJI{i}(j).s.u.p.c      = objI{i}(j).s.u.p.c;
                            if ~isempty(objI{i}(j).s.u.r)
                                fprintf('Cycle #%3d\n',j-1);
                            end
                        else
                            break;
                        end
                    end
                end
                if size(inp.C,2) ~= 1
                    save(SD,"OBJI");
                end
                OBJ2 = OBJI;
            else
                VD   = dir(SD);
                OBJ1 = struct2cell(load(VD.name));
                OBJ2 = cat(2,OBJ1{:});
            end
            %  > Plot...
            switch size(inp.C,2)
                case 1
                    for i = 1:numel(OBJ2)
                        PP_1D.Plot(inp,msh,OBJ2{i}(end),...
                            {[0,0,0,0,0],...
                             [0,0],...
                             [0],...
                             [0,0],...
                             [v(i)-1 == 1,v(i)-1 == 2],...
                             [0],...
                             [0],...
                             [0,0,0,0]},{});
                    end
                    PP_1D.Plot(inp,msh,OBJ2,...
                            {[0,0,0,0,0],...
                             [0,0],...
                             [0],...
                             [0,0],...
                             [0,0],...
                             [1],...
                             [0],...
                             [0,0,0,0]},{v});
                case 2
                    % for i = 1:numel(OBJ2)
                    %     PP_2D.Plot(inp,msh,OBJ2{i}(end),...
                    %         {[0],...
                    %          [0,0],...
                    %          [0,0],...
                    %          [1],...
                    %          [0,0],...
                    %          [0,0,0]});
                    % end
            end
            P_Ad.Plot(inp,msh,OBJ2,v);
        end
        % >> 1.2. ---------------------------------------------------------
        %  > Check stopping criteria.
        function [f] = Check_SC(inp,obj,u)
            switch u
                case 1, rp = obj(end).e.a   .n_abs.c(2,1);
                case 2, rp = obj(end).e.x(2).n_abs.c(2,1);
                case 3, rp = obj(end).e.x(1).n_abs.c(2,1);
            end
            c = [isempty(obj(end).s.u.s.c{1}),...
                 inp.P_Ad.ec   >= rp,...
                 inp.P_Ad.fNNZ <=  obj(end).s.c.m.nnz.At./obj(1).s.c.m.nnz.At,...
                 inp.P_Ad.nc   <= size(obj,2)-1];
             if c(1), fprintf('Stopping criterion: no faces left to be added.\n');
             end
             if c(2), fprintf('Stopping criterion: min. error treshold (L_inf norm).\n');
             end
             if c(3), fprintf('Stopping criterion: max. number of nonzeros (or DoFs).\n');
             end
             if c(4), fprintf('Stopping criterion: max. number of cycles.\n');
             end
             f = any(c,2);
        end
        
        %% > 2. -----------------------------------------------------------
        % >> 2.1. ---------------------------------------------------------
        %  > Update fields 'e' and 's'.
        function [obj] = Update_P(inp,msh,e,f,s,u)
            %  > Auxiliary variables.
            fn = reshape(setdiff(string(fieldnames(s)),"u"),1,[]);
            fD = true   (size(fn));
            rl = unique (ceil(s.u.p.c{1}./2)','rows')';
            ra = true   (size(rl));
            ss = s;
            
            while 1
                %  > Select face(s)...
                switch u
                    case 1, rp = e.a   .t.f_abs(:,3);
                    case 2, rp = e.x(2).t.f_abs(:,3);
                    case 3, rp = e.x(1).t.f_abs(:,1);
                    otherwise
                        return;
                end
                F = P_Ad.fCR(inp,round(rp,10),ra);
                if ~isempty(F.r)
                    %  > Refinement rules.
                    fA = cell(size(F.r));
                    for i = 1:size(F.r,1)
                        fA{i,1} = P_Ad.RR(inp,F.r(i,1),s.c.s.sp.f,s.c.s.ss.f,rl);
                    end
                    fB = RunLength(sort(cat(1,fA{:,1}),1));
                    %  > Mesh regularisation (w/ rule verification).
                    fC = P_Ad.MR(inp,msh,fB,s.c.s.sp,rl);
                    
                    %  > Update field 's'.
                    for i = fn
                        for j = 1:numel(s.u.p.(i))
                            s.u.p.(i){j}(fC,:) = 2.*ceil(s.u.p.(i){j}(fC,:)./2)+1;
                            s.u.s.(i){j}       = fC;
                        end
                        if i == "c"
                            s.u.r = F.r;
                        end
                    end
                    switch size(inp.C,2)
                        case 1, s = B1_1D.Update_s(inp,msh,f,s);
                        case 2, s = B1_2D.Update_s(inp,msh,f,s);
                    end
                    for i = 1:numel(fn)
                        fD(i) = ~any(any(cellfun(@isempty,s.(fn(i)).s.M),1),2);
                    end
                    if all(fD,2)
                        %  > Update field 'e'.
                        e = NE_XD.Update_e(msh,e,s);
                        break;
                    else
                        s = ss;
                        for i = fn
                            for j = 1:numel(s.u.p.(i))
                                s.u.s.(i){j} = [];
                            end
                            if i == "c"
                                s.u.r = [];
                            end
                        end
                        break;
                    end
                else
                    for i = fn
                        for j = 1:numel(s.u.p.(i))
                            s.u.s.(i){j} = F.r;
                        end
                    end
                    break;
                end
            end
            %  > Assign to 'obj'...
            obj.e = e;
            obj.f = f;
            obj.s = s;
        end
        % >> 2.2. ---------------------------------------------------------
        %  > Flag faces to be coarsened/refined.
        function [F] = fCR(inp,e,ra)
            % >> Coarsening.
            %  > \NA.
            % >> Refinement.
            F.r = P_Ad.fR(inp,e,ra);
        end
        %  > 2.2.1. -------------------------------------------------------
        %  > Auxiliary function #1 (coarsening).
        %  > \NA.
        %  > 2.2.2. -------------------------------------------------------
        %  > Auxiliary function #2 (refinement).
        function [A] = fR(inp,e,ra)
            %  > Select faces (fixed-fraction parameter approach).
            A = [];
            if any(ra,1)
                ec = sort(e(ra),1,"descend");
                A  = find(e >= ec(ceil(inp.P_Ad.fp.*sum(ra,1)),1) & ra);
            end
        end
        
        %% > 3. -----------------------------------------------------------
        % >> Isotropic coarsening routines.
        %  > \NA.
        
        %% > 4. -----------------------------------------------------------
        % >> Isotropic refinement routines.
        % >> 4.1. ---------------------------------------------------------
        %  > 4.1.1. -------------------------------------------------------
        %  > Refinement rules.
        function [D] = RR(inp,A,np,ns,rl)
            % >> Rule #1 (neighbouring rule).
            if any(rl(np{A,1},1) < rl(A,1),1)
                B = np{A,1}(rl(np{A,1},1)  < rl(A,1));
            else
                B = np{A,1}(rl(np{A,1},1) == rl(A,1));
            end
            %  > Check whether faces in 'B' are eligible for refinement. If not, refine the ones on the secondary layer/zone.
            C = P_Ad.loop(B,np,rl);
            % >> Rule #2 (irregular rule).
            if ~inp.P_Ad.R2
                D = C;
            else
                D = P_Ad.loop(C,ns,rl);
            end
        end
        %  > 4.1.2. -------------------------------------------------------
        %  > Auxiliary function #1.
        function [E] = loop(A,n,rl)
            %  > Initialize...
            D = cell (size(A));
            k = false(size(A));
            
            for i = 1:size(k,1)
                %  > If any(X,1), analyse face.
                %    Else, move on.
                k(i,1) = any(rl(n{A(i,1),1},1) < rl(A(i,1),1),1);
                if k(i,1)
                    %  > Neighbouring elements (or faces) in 'A(i,1)' w/ rl<rl*.
                    B = n{A(i,1),1}(rl(n{A(i,1),1},1) < rl(A(i,1),1),1);
                    %  > Refine neighbouring elements (or faces) until the ones we wish to refine are eligible for refinement.
                    while 1
                        C = cell (size(B));
                        l = false(size(B));
                        for j = 1:size(B,1)
                            if any(rl(n{B(j,1),1},1) < rl(B(j,1),1),1)
                                C{j,1} = n{B(j,1),1}(rl(n{B(j,1),1},1) < rl(B(j,1),1),1);
                                l(j,1) = true;
                            end
                        end
                        %  > If all elements (or faces) in 'B' are eligible for refinement, refine them.
                        %    Else, re-do loop until its lower-order neighbouring elements (or faces) have been refined.
                        if all(~l,1)
                            break;
                        else
                            B = reshape(unique(cat(1,C{:,1})),[],1);
                        end
                    end
                    D{i,1} = B;
                end
            end
            if ~any(k,1)
                E = A;
            else
                E = reshape(unique(cat(1,D{:,1})),[],1);
            end
        end
        % >> 4.2. ---------------------------------------------------------
        % >> Mesh regularisation (extend refinement to avoid the creation of lower-order regions in between higher-order ones).
       function [B] = MR(inp,msh,A,n,rl)
            %  > Initialize...
            B  = A;
            rv = rl;
            
            while 1
                %  > Flag neighbouring elements (or faces).
                C       = AuxiliaryFunctions.setdiff(RunLength(sort(cat(1,n.f{A,1}),1)),B);
                %       = AuxiliaryFunctions.setdiff(RunLength(sort(reshape(msh.c.f.if(reshape(cat(2,msh.f.ic{A,1}),[],1),:),[],1),1)),B);
                %  > Update (virtual) level of refinement.
                rv(A,1) = rv(A,1)+1;
                rt      = rv;
                %  > Tag neighbouring elements (or faces)...
                j = false(size(C));
                for i = 1:size(C,1)
                    j(i,1) = all(any(reshape(rv(msh.c.f.if(n.c{C(i,1),1},:),1),size(n.c{C(i,1),1},1),[]) > rv(C(i,1),1),2) | all(ismembc(msh.c.f.if(n.c{C(i,1),1},:),C),2),1);
                end
                if any(j,1)
                    %  > ... and decide whether to refine them or not.
                    D       = C (j,1);
                    rt(D,1) = rt(D,1)+1;
                    k       = false(size(D));
                    for i = 1:size(D,1)
                        k(i,1) = any(all(reshape(rt(msh.c.f.if(n.c{D(i,1),1},:),1),size(n.c{D(i,1),1},1),[]) == rt(D(i,1),1),2),1);
                    end
                    if all(~k,1)
                        break;
                    else
                        [A,B] = P_Ad.VR(inp,msh,B,D(k,1),rv,n);
                    end
                else
                    break;
                end
            end
        end
        %  > 4.2.1. -------------------------------------------------------
        %  > Auxiliary function #1.
        function [A,B] = VR(inp,msh,B,D,rv,np)
            %  > Initialize...
            j       = false(size(D));
            rv(D,1) = rv(D,1)+1;
            
            if ~inp.P_Ad.R2
                %  > Verify rule  #1 (only).
                for i = 1:size(D,1)
                    j(i,1) = all(range(reshape(rv(msh.c.f.if(np.c{D(i,1),1},:),1),size(np.c{D(i,1),1},1),[]),2) <= 1,1);
                end
            else
                %  > Verify rules #1 and #2.
                for i = 1:size(D,1)
                    if max(rv(np.f{D(i,1),1},1))-min(rv(np.f{D(i,1),1},1)) <= 1
                        j(i,1) = true;
                    end
                end
            end
            D = D(j,1);
            %  > Assign...
            B = sort(cat(1,B,D),1);
            A = D;
        end
        
        %% > 5. -----------------------------------------------------------
        % >> 5.1. ---------------------------------------------------------
        %  > Plot(s).
        function [] = Plot(inp,msh,obj,v)
            % >> ----------------------------------------------------------
            % >> Auxiliary variables.
            %  > ----------------------------------------------------------
            F            = [1,1];
            U            = [1];
            objP.C       = [linspecer(9,'Qualitative')]; %  > Color.
            objP.C(10,:) = [0,0,0];                      %  > k.
            objP.C(11,:) = [1,1,1];                      %  > w.
            objP.E       = [0];                          %  > Export?
            objP.MS      = [8.50];                       %  > Marker size.
            objP.S       = ["1","\infty"];               %  > Lk error norm(s).
            if objP.E
                %LF = zeros(numel(F),2);
                for i = 1:numel(F)
                    for j = [1,2]
                        LF(i,j) = join(["AMR (1.",num2str(i),".",num2str(j),").pdf"],'');
                    end
                end
                switch size(inp.C,2)
                    case 1, PF = "[Post-processing]/[.PDF Files]/[1D]";
                    case 2, PF = "[Post-processing]/[.PDF Files]/[2D]";
                end
            end
            %  > ----------------------------------------------------------
            L                         = {'$\Phi$',...
                                         '$\phi_{\left(h\right)}$',...
                                         '$\phi_{\left(\hspace{0.085em}l\hspace{0.085em}\right)}$'};
            M                         = [":o","-.v","-.^"];
            objP.Axes.Clipping        = [1];
            objP.Axes.FontSize{1}     = [25.00,25.00];
            objP.Axes.FontSize{2}     = [32.50,32.50];
            switch size(inp.C,2)
                case 1, objP.Axes.Label(1) = "$\textrm{Number of DoFs (ratio)}$";
                case 2, objP.Axes.Label(1) = "$\sqrt{\textrm{Number of DoFs (ratio)}}$";
                otherwise
                    return;
            end
            objP.Axes.Reverse         = [0,0];
            objP.Axes.Scale           = [1,1];
            objP.Axes.Tick.X          = [0];
            objP.Axes.Tick.XTick      = [1,2,3,4,5,10];
            objP.Axes.Tick.XTickLabel = [];
            objP.Axes.Tick.Y          = [0];
            objP.Axes.Tick.YTick      = [];
            objP.Axes.Tick.YTickLabel = [];
            objP.ET                   = [1e-10];
            objP.Grid.Display         = [1];
            objP.Legend.Display       = [1];
            objP.Legend.FontSize      = [30.00];
            objP.Legend.Location      = ['NorthEast'];
            objP.Legend.NumColumns    = [1];
            objP.LW                   = [3.50,1.00];
            objP.pbaspect             = [1,1,1];
            %  > ----------------------------------------------------------
            [nnz,NNZ,Y] = P_Ad.AF_1(inp,obj,v);
            if U
                [NNZU,YU] = P_Ad.AF_2(inp,msh,nnz);
            end
            %  > ----------------------------------------------------------
            if U
                A  = repmat(NNZU(:,1),1,2);
                B  = cell(numel(F),2);
                pu = cell(numel(F),2);
                for i = 1:numel(F)
                    for j = [1,2]
                        for k = 1:size(NNZU,1)
                            for l = 1:size(NNZU,2)
                                pu{i,j}(k,l) = YU{i,l}(j,k);
                            end
                            B{i,j}(k,:) = [YU{i,1}(j,k),objP.ET];
                        end
                    end
                end
                if ~objP.E
                    if size(NNZU,1) <= 5
                        PFLex = [-105,-62.5];
                    else
                        PFLex = [-120,-67.5];
                    end
                else
                    if size(NNZU,1) <= 5
                        PFLex = [-407.5,-10];
                    else
                        PFLex = [-422.5,-10];
                    end
                end
                for i = 1:size(NNZU,1)
                    objP.Legend.Label(i)  = [join(["$\mathbf{\hat{p}}_{",num2str(2.*i-1),"}$"],'')];
                end
            else
                objP.Legend.Label = L(v);
            end
            %  > ----------------------------------------------------------
            PU = cell(1,size(NNZU,1));
            Z  = cell(1,size(Y,1));
            for i = 1:numel(F)
                if F(i)
                    for j = [1,2]
                        %  > ----------------------------------------------
                        switch i
                            case 1, objP.Axes.Label(2) = [join(["$\|\bar{\tau}_{f_{\left(k\right)}}\|_{_{",objP.S(j),"}}$"],'')];
                            case 2, objP.Axes.Label(2) = [join(["$\|e_{c}^{\hspace{0.085em}}\|_{_{",objP.S(j),"}}$"],'')];
                        end
                        for k = 1:size(Y,1)
                            Z{k} = Y{k,i}(j,:);
                        end
                        objP.Axes.Limits = [PP_1D.AF_2(objP,NNZ,[],Z,[])];
                        %  > ----------------------------------------------
                        if ~objP.E
                            if j == 1
                                fprintf('Plotting...\n');
                                figure ('Color',objP.C(11,:),'Windowstate','Maximized');
                            end
                            subplot(1,2,j);
                            hold on;
                        else
                            fprintf('Plotting...\n');
                            figure ('Color',objP.C(11,:),'Windowstate','Maximized');
                            hold on;
                        end
                        if ~U
                            for k = size(NNZ,1):-1:1
                                PA{k} = plot(NNZ{k},Y{k,i}(j,:),M(v(k)),'Color',objP.C(10,:),'LineWidth',objP.LW(2),'MarkerFaceColor',objP.C(11,:),'MarkerSize',objP.MS(1));
                            end
                            PP_1D.AF_1(objP,PA);
                        else
                            for k = 1:size(NNZU,1)
                                plot(A(k,:),B{i,j}(k,:),"--",'Color',objP.C(k,:),'LineWidth',objP.LW(2));
                                if size(NNZU,2) == 1
                                    PU{k} = plot(NNZU(k,:),pu{i,j}(k,:), "o",'Color',objP.C(k,:),'MarkerFaceColor',objP.C(k,:),'LineWidth',objP.LW(1));
                                else
                                    PU{k} = plot(NNZU(k,:),pu{i,j}(k,:),"-o",'Color',objP.C(k,:),'MarkerFaceColor',objP.C(k,:),'LineWidth',objP.LW(1));
                                end
                            end
                            for k = size(NNZ,1):-1:1
                                PA{k} = plot(NNZ{k},Y{k,i}(j,:),M(v(k)),'Color',objP.C(10,:),'LineWidth',objP.LW(2),'MarkerFaceColor',objP.C(11,:),'MarkerSize',objP.MS(1));
                            end
                            legendflex(cat(2,PA{:}),L(v),'Anchor',{'ne','ne'},'Buffer',PFLex,'Fontsize',objP.Legend.FontSize);
                            PP_1D.AF_1(objP,PU);
                        end
                        if objP.E
                            exportgraphics(gcf,LF(i,j),'ContentType','Vector');
                            movefile      (    LF(i,j),PF);
                        end
                        %  > ----------------------------------------------
                    end
                end
            end
            %  > ----------------------------------------------------------
        end
        %  > 5.1.1. -------------------------------------------------------
        %  > Auxiliary function #1.
        function [nnz,NNZ,Y] = AF_1(inp,obj,v)
            for i = 1:numel(obj)
                switch v(i)
                    case 1
                        Y(i,:) = ...
                            {[arrayfun(@(x) x.e.a   .n_abs.t.f(1,3),obj{i});
                              arrayfun(@(x) x.e.a   .n_abs.t.f(2,3),obj{i})];
                             [arrayfun(@(x) x.e.a   .n_abs.c  (1,1),obj{i});
                              arrayfun(@(x) x.e.a   .n_abs.c  (2,1),obj{i})]};
                    case 2
                        Y(i,:) = ...
                            {[arrayfun(@(x) x.e.x(2).n_abs.t.f(1,3),obj{i});
                              arrayfun(@(x) x.e.x(2).n_abs.t.f(2,3),obj{i})];
                             [arrayfun(@(x) x.e.x(2).n_abs.c  (1,1),obj{i});
                              arrayfun(@(x) x.e.x(2).n_abs.c  (2,1),obj{i})]};
                    case 3
                        Y(i,:) = ...
                            {[arrayfun(@(x) x.e.x(1).n_abs.t.f(1,1),obj{i});
                              arrayfun(@(x) x.e.x(1).n_abs.t.f(2,1),obj{i})];
                             [arrayfun(@(x) x.e.x(1).n_abs.c  (1,1),obj{i});
                              arrayfun(@(x) x.e.x(1).n_abs.c  (2,1),obj{i})]};
                    otherwise
                        return;
                end
                nnz{i,1} = arrayfun(@(x) x.s.c.m.nnz.At,obj{i}).^(1./size(inp.C,2));
                NNZ{i,1} = nnz{i,1}./nnz{i,1}(1);
            end
        end
        %  > 5.1.2. -------------------------------------------------------
        %  > Auxiliary function #2.
        function [NNZU,YU] = AF_2(inp,msh,nnz)
            %  > Auxiliary variables.
            if size(inp.C,2) == 1
                LD(1) = 1;
                LD(2) = 1;
            else
                LD(1) = 1;
                LD(2) = 1;
                VS    = "[P-Adaptation]/[.mat Files]/[2D]/[.T2]/A2_2D.mat";
            end
            
            if ~LD(1)
                m = 1;
                n = 2;
                for i = 1:m
                    mshU(i) = msh; %#ok<*AGROW>
                    for j = 1:n
                        if i == 1
                            inp.PS.EE    = 0;
                            inpU(j)      = inp;
                            inpU(j).PS.P = 2.*(6+j)-1;
                        end
                        switch size(inp.C,2)
                            case 1
                                objU(i,j).f = A3_1D.Initialize_f(inpU(j),mshU(i));
                                objU(i,j).s = B1_1D.Initialize_s(inpU(j),mshU(i));
                                objU(i,j).s = B1_1D.Update_s    (inpU(j),mshU(i),objU(i,j).f,objU(i,j).s);
                            case 2
                                objU(i,j).f = A3_2D.Initialize_f(inpU(j),mshU(i));
                                objU(i,j).s = B1_2D.Initialize_s(inpU(j),mshU(i));
                                objU(i,j).s = B1_2D.Update_s    (inpU(j),mshU(i),objU(i,j).f,objU(i,j).s);
                        end
                        objU(i,j).e            = NE_XD.Initialize_e(inpU(j),mshU(i));
                        objU(i,j).e            = NE_XD.Update_e    (        mshU(i),objU(i,j).e,objU(i,j).s);
                        OBJU(i,j).e            = objU(i,j).e;
                        OBJU(i,j).s.c.m.nnz.At = objU(i,j).s.c.m.nnz.At;
                        fprintf('#%2d\n',j);
                    end
                end
                if size(inp.C,2) ~= 1
                    save(VS,"OBJU");
                end
                OBJ2 = OBJU;
            else
                D    = dir(VS);
                OBJ1 = struct2cell(load(D.name));
                OBJ2 = cat(2,OBJ1{:}); 
            end
            for i = 1:size(OBJ2,1)
                YU(:,i) = ...
                    {[arrayfun(@(x) x.e.a.n_abs.t.f(1,3),OBJ2(i,:));
                      arrayfun(@(x) x.e.a.n_abs.t.f(2,3),OBJ2(i,:))];
                     [arrayfun(@(x) x.e.a.n_abs.c  (1,1),OBJ2(i,:));
                      arrayfun(@(x) x.e.a.n_abs.c  (2,1),OBJ2(i,:))]};
            end
            NNZU = (arrayfun(@(x) x.s.c.m.nnz.At,OBJ2).^(1./size(inp.C,2)))'./nnz{1}(1);
            if LD(2)
                NNZUM = 15;
                NNZU  = cat(2,NNZU,repmat(NNZUM,size(NNZU,1),1));
                for i = 1:size(YU,1)
                    if ~(i == 1 && size(inp.C,2) == 2)
                        p = 2;
                    else
                        p = 3;
                    end
                    YU{i,size(NNZU,2)} = YU{i,size(NNZU,2)-1}.*exp(-p.*(1:size(NNZU,1)).*log(NNZU(:,end)./NNZU(:,end-1))');
                end
            end
        end
    end
end