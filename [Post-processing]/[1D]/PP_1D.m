classdef PP_1D
    methods (Static)
        %% > 1. -----------------------------------------------------------
        % >> 1.1. ---------------------------------------------------------
        %  > 1D plot(s).
        function [] = Plot(inp,msh,obj,F,Z)
            % >> ----------------------------------------------------------
            % >> Auxiliary variables.
            %  > ----------------------------------------------------------
            ms           = [1];                          %  > Manufactured solution (MMS).
            objP.C       = [linspecer(9,'Qualitative')]; %  > Color.
            objP.C(10,:) = [0,0,0];                      %  > k.
            objP.C(11,:) = [1,1,1];                      %  > w.
            objP.E       = [0];                          %  > Export?
            objP.MS      = [6.50];                       %  > Marker size.
            objP.S       = ["1","\infty"];               %  > Lk error norm(s).
            if size(msh,2) == 1
                mX = min(msh.f.Xv(:,1),[],1);
                MX = max(msh.f.Xv(:,1),[],1);
            end
            if ~inp.PS.EE
                F{1}([2:5]) = 0;
                F{4}([1:2]) = 0;
                F{5}([1:2]) = 0;
            end
            if all(cellfun(@isempty,Z),2)
                F{2}([1:2]) = 0;
                F{3}([1])   = 0;
            end
            if size(obj,1) == 1
                F{8}([1:4]) = 0;
            else
                F{1}([1:5]) = 0;
                F{2}([1:2]) = 0;
                F{3}([1])   = 0;
                F{4}([1:2]) = 0;
                F{5}([1:2]) = 0;
                F{6}([1])   = 0;
                F{7}([1])   = 0;
            end
            %  > ----------------------------------------------------------
            if objP.E
                LF = cell(numel(F),1);
                for i = 1:numel(F)
                    switch i
                        case 1
                            for j = 1:numel(F{i})
                                LF{i}(j) = join(["NEE (1.",num2str(j-1),").pdf"],'');
                            end
                        case 2
                            for j = 1:numel(F{i})
                                LF{i}(j) = join(["AD (1.1.",num2str(j),").pdf"],'');
                                %         = join(["AD (1.2.",num2str(j),").pdf"],'');
                                %         = join(["AD (1.3.",num2str(j),").pdf"],'');
                            end
                        case 3
                            for j = 1:numel(Z)
                                LF{i}(j) = join(["AD (2.",num2str(j),").pdf"],'');
                            end
                        case 4
                            for j = 1:numel(F{i})
                                LF{i}(j) = join(["AD (3.",num2str(j),").pdf"],'');
                            end
                        case 5
                            for j = 1:numel(F{i})
                                LF{i}(j) = join(["AMR (1.3.",num2str(j),").pdf"],'');
                            end
                        case 6
                            LF{i} = "AMR (1.4).pdf";
                        case 7
                            for j = [1,2]
                                LF{i}(j) = join(["B (1.",num2str(j),").pdf"],'');
                            end
                        case 8
                            switch ms
                                case 1
                                    A = [1,1;1,1;3,3;3,3];
                                    B = [1,2;3,4;1,2;3,4];
                                case 2
                                    A = [2,2;2,2;4,4;4,4];
                                    B = [1,2;3,4;1,2;3,4];
                            end
                            for j = 1:numel(F{i})
                                for k = [1,2]
                                    LF{i}(j,k) = join(["NEE (D",num2str(A(j,k)),".",num2str(B(j,k)),").pdf"],'');
                                end
                            end
                        otherwise
                            return;
                    end
                end
            end
            PF = "[Post-processing]/[.PDF Files]/[1D]";
            %  > ----------------------------------------------------------
            
            % >> ----------------------------------------------------------
            % >> Plot...
            %  > ----------------------------------------------------------
            if any(F{1},2)
                %  > ------------------------------------------------------
                objP.Axes.Clipping        = [1];
                objP.Axes.FontSize{1}     = [25.00,25.00];
                objP.Axes.FontSize{2}     = [37.50,32.50];
                objP.Axes.Label           = ["$x$","$\textrm{Error magnitude}$"];
                objP.Axes.Reverse         = [0,0];
                objP.Axes.Scale           = [0,1];
                objP.Axes.Tick.X          = [0];
                objP.Axes.Tick.XTick      = [linspace(mX,MX,5)];
                objP.Axes.Tick.XTickLabel = [];
                objP.Axes.Tick.Y          = [0];
                objP.Axes.Tick.YTick      = [];
                objP.Axes.Tick.YTickLabel = [];
                objP.ET                   = [1e-10];
                objP.Grid.Display         = [1];
                objP.Legend.Display       = [1];
                objP.Legend.FontSize      = [30.00];
                objP.Legend.Location      = ['NorthEast'];
                objP.LW                   = [3.50,2.50];
                objP.pbaspect             = [1,1,1];
                L                         = [cell(1,numel(F{1}))];
                M                         = [cell(1,numel(F{1}))];
                N                         = [cell(1,numel(F{1}))];
                X                         = [cell(1,numel(F{1}))];
                Y                         = [cell(1,numel(F{1}))];
                %  > ------------------------------------------------------
                for i = 1:numel(F{1})
                    if F{1}(i)
                        %  > ----------------------------------------------
                        switch i
                            case 1
                                L{i} = ["$|\bar{\tau}_{f}^{\Phi}|$",...
                                        "$|\bar{\tau}_{c}^{\Phi}|$",...
                                        "$|e_{c}^{\hspace{0.085em}\Phi}|$"];
                                M{i} = ["-o","-o","-o"];
                                N{i} = [1,2,3];
                                X{i} = {msh.f.Xv,msh.c.Xc,msh.c.Xc};
                                Y{i} = {obj(end).e.a.t.f_abs(:,1),obj(end).e.a.t.c_abs(:,1),obj(end).e.a.c.c_abs(:,1)};
                            case 2
                                L{i} = ["$|\bar{\tau}_{f}^{\phi_{\left(h\right)}}|$",...
                                        "$|\bar{\tau}_{f_{\left(0\right)}}^{\phi_{\left(h\right)}}|$",...
                                        "$|\bar{\tau}_{f_{\left(1\right)}}^{\phi_{\left(h\right)}}|$",...
                                        "$|\bar{\tau}_{f}^{\Phi-\phi_{\left(h\right)}}|$",...
                                        "$|\bar{\tau}_{f_{\left(0\right)}}^{\Phi-\phi_{\left(h\right)}}|$",...
                                        "$|\bar{\tau}_{f_{\left(1\right)}}^{\Phi-\phi_{\left(h\right)}}+\ldots|$"];
                                M{i} = [":o","--o","--o";":","--","--"];
                                N{i} = [1,4,5;1,4,5];
                                X{i} = {msh.f.Xv,msh.f.Xv,msh.f.Xv;
                                        msh.f.Xv,msh.f.Xv,msh.f.Xv};
                                Y{i} = {obj(end).e.x(2).t.f_abs(:,1),obj(end).e.x(2).t.f_abs(:,2),obj(end).e.x(2).t.f_abs(:,3);
                                        obj(end).e.d(2).t.f_abs(:,1),obj(end).e.d(2).t.f_abs(:,2),obj(end).e.d(2).t.f_abs(:,3)};
                            case 3
                                L{i} = ["$|\bar{\tau}_{c}^{\phi_{\left(h\right)}}|$",...
                                        "$|e_{c}^{\hspace{0.085em}\phi_{\left(h\right)}}|$",...
                                        "$|\bar{\tau}_{c}^{\Phi-\phi_{\left(h\right)}}|$",...
                                        "$|e_{c}^{\hspace{0.085em}\Phi-\phi_{\left(h\right)}}|$"];
                                M{i} = ["-o","-o";"-","-"];
                                N{i} = [2,3;2,3];
                                X{i} = {msh.c.Xc,msh.c.Xc;
                                        msh.c.Xc,msh.c.Xc};
                                Y{i} = {obj(end).e.x(2).t.c_abs(:,1),obj(end).e.x(2).c.c_abs(:,1);
                                        obj(end).e.d(2).t.c_abs(:,1),obj(end).e.d(2).c.c_abs(:,1)};
                            case 4
                                L{i} = ["$|\bar{\tau}_{f}^{\phi_{\left(l\right)}}|$",...
                                        "$|\bar{\tau}_{f_{\left(0\right)}}^{\Phi}\hspace{0.095em}|$",...
                                        "$|\bar{\tau}_{f_{\left(1\right)}}^{\phi_{\left(l\right)}}|$",...
                                        "$|\bar{\tau}_{f}^{\Phi-\phi_{\left(l\right)}}|$",...
                                        "",...
                                        "$|\bar{\tau}_{f_{\left(1\right)}}^{\Phi-\phi_{\left(l\right)}}+\ldots|$"];
                                M{i} = [":o","--o","--o";":","--","--"];
                                N{i} = [1,4,5;1,11,5];
                                X{i} = {msh.f.Xv,msh.f.Xv,msh.f.Xv;
                                        msh.f.Xv,msh.f.Xv,msh.f.Xv};
                                Y{i} = {obj(end).e.x(1).t.f_abs(:,1),obj(end).e.a.t.f_abs(:,2),obj(end).e.x(1).t.f_abs(:,1);
                                        obj(end).e.d(1).t.f_abs(:,1),NaN(msh.f.Nf,1)          ,obj(end).e.d(1).t.f_abs(:,3)};
                            case 5
                                L{i} = ["$|\bar{\tau}_{c}^{\phi_{\left(l\right)}}|$",...
                                        "$|e_{c}^{\hspace{0.085em}\phi_{\left(l\right)}}|$",...
                                        "$|\bar{\tau}_{c}^{\Phi-\phi_{\left(l\right)}}|$",...
                                        "$|e_{c}^{\hspace{0.085em}\Phi-\phi_{\left(l\right)}}|$"];
                                M{i} = ["-o","-o";"-","-"];
                                N{i} = [2,3;2,3];
                                X{i} = {msh.c.Xc,msh.c.Xc;
                                        msh.c.Xc,msh.c.Xc};
                                Y{i} = {obj(end).e.x(1).t.c_abs(:,1),obj(end).e.x(1).c.c_abs(:,1);
                                        obj(end).e.d(1).t.c_abs(:,1),obj(end).e.d(1).c.c_abs(:,1)};
                        end
                        %  > ----------------------------------------------
                        objP.Axes.Limits  = [PP_1D.AF_2(objP,X{i},[],Y{i},[])];
                        objP.Legend.Label = [L{i}];
                        if mod(numel(Y{i}),2)
                            objP.Legend.NumColumns = [1];
                        else
                            objP.Legend.NumColumns = [2];
                        end
                        %  > ----------------------------------------------
                        fprintf('Plotting...\n');
                        figure ('Color',objP.C(11,:),'Windowstate','Maximized');
                        hold on;
                        PP_1D.P1(...
                            [objP],...
                            [M{i}],...
                            [N{i}],...
                            [X{i}],...
                            [Y{i}]);
                        if objP.E
                            exportgraphics(gcf,LF{1}(i),'ContentType','Vector');
                            movefile      (    LF{1}(i),PF);
                        end
                        %  > ----------------------------------------------
                    end
                end
                %  > ------------------------------------------------------
            end
            %  > ----------------------------------------------------------
            %  > #2.
            if F{2}
                %  > ------------------------------------------------------
                objP.Axes.Clipping        = [1];
                objP.Axes.FontSize{1}     = [25.00,25.00];
                objP.Axes.FontSize{2}     = [37.50,32.50];
                objP.Axes.Label           = ["$x$","$\textrm{Error magnitude}$"];
                objP.Axes.Reverse         = [0,0];
                objP.Axes.Scale           = [0,1];
                objP.Axes.Tick.X          = [0];
                objP.Axes.Tick.XTick      = [linspace(mX,MX,5)];
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
                objP.LW                   = [3.50];
                objP.pbaspect             = [1,1,1];
                x                         = [4.337906412938963e-04,...
                                             1.279699999998663e-04]; %#ok<NASGU>
                L                         = [cell(1,numel(F{2}))];
                M                         = [cell(1,numel(F{2}))];
                N                         = [cell(1,numel(F{2}))];
                X                         = [cell(1,numel(F{2}))];
                Y                         = [cell(1,numel(F{2}))];
                %  > ------------------------------------------------------
                for i = 1:numel(F{2})
                    if F{2}(i)
                        %  > ----------------------------------------------
                        switch i
                            case 1
                                L{i} = ["$|\bar{\tau}_{f}^{\Phi}\hspace{0.800em}|$",...
                                        "$|\bar{\tau}_{f_{\left(\Sigma k\right)}}^{\Phi}|$"];
                                M{i} = ["-o","-o"];
                                N{i} = [1,5;1,5];
                                X{i} = {msh.f.Xv,msh.f.Xv};
                                Y{i} = {obj.e.a.t.f_abs(:,1),obj.e.a.t.f_abs(:,3)};
                            case 2
                                L{i} = ["$|\bar{\tau}_{c}^{\Phi}|$",...
                                        "$|e_{c}^{\hspace{0.085em}\Phi}|$"];
                                M{i} = ["-o","-o"];
                                N{i} = [1,5;2,3];
                                X{i} = {msh.c.Xc,msh.c.Xc};
                                Y{i} = {obj.e.a.t.c_abs(:,1),obj.e.a.c.c_abs(:,1)};
                        end
                        %  > ----------------------------------------------
                        objP.Axes.Limits  = [PP_1D.AF_2(objP,X{i},[],Y{i},[])];
                        objP.Legend.Label = [L{i}];
                        %  > ----------------------------------------------
                        fprintf('Plotting...\n');
                        figure ('Color',objP.C(11,:),'Windowstate','Maximized');
                        hold on;
                        PP_1D.P2(...
                            [objP],...
                            [M{i}],...
                            [N{i}],...
                            [X{i}],...
                            [Y{i}],[Z]);
                        %  if i == 1
                        %      plot(msh.f.Xv,repmat(x(1),size(msh.f.Xv,1),1),"--",'Color',objP.C(1,:),'LineWidth',1.00);
                        %      plot(msh.f.Xv,repmat(x(2),size(msh.f.Xv,1),1),"--",'Color',objP.C(5,:),'LineWidth',1.00);
                        %  end
                        if objP.E
                            exportgraphics(gcf,LF{2}(i),'ContentType','Vector');
                            movefile      (    LF{2}(i),PF);
                        end
                        %  > ----------------------------------------------
                    end
                end
                %  > ------------------------------------------------------
            end
            %  > ----------------------------------------------------------
            %  > #3.
            if F{3}
                %  > ------------------------------------------------------
                Yb = [round(max(cat(2,Z{:}),[],'all'),1)];
                if max(cat(2,Z{:}),[],'all') > Yb
                    Yb = Yb+0.1;
                end
                Yc = [RunLength(cat(2,0:0.1:Yb,1))];
                for i = 1:numel(Yc)
                    objP.Axes.Tick.YTick(i) = Yc(i);
                    objP.Axes.Tick.YTickLabel{i} = join(["$",num2str(objP.Axes.Tick.YTick(i)),"$"],'');
                    if numel(num2str(objP.Axes.Tick.YTick(i))) == 1
                        objP.Axes.Tick.YTickLabel{i} = join(["$",num2str(objP.Axes.Tick.YTick(i)),".0$"],'');
                    end
                end
                if objP.Axes.Tick.YTick(end-1)+0.1 < 1-0.1
                    objP.Axes.Tick.YTick(end) = objP.Axes.Tick.YTick(end-1)+2.*0.1;
                end
                %
                %  > ------------------------------------------------------
                objP.Axes.Clipping        = [0];
                objP.Axes.FontSize{1}     = [23.50,25.00];
                objP.Axes.FontSize{2}     = [32.50,32.50];
                objP.Axes.Limits.X        = [ 0.50, 2.50];
                objP.Axes.Limits.Y        = [min(objP.Axes.Tick.YTick,[],2),max(objP.Axes.Tick.YTick,[],2)];
                objP.Axes.Reverse         = [0,0];
                objP.Axes.Scale           = [0,0];
                objP.Axes.Tick.XTick      = [1,2];
                objP.Axes.Tick.XTickLabel = ["$\textrm{Criterion:}|\bar{\tau}_{f}^{\Phi}|\textrm{ (tau)}$";...
                                             "$\textrm{Criterion:}|\bar{\tau}_{f_{\left(\Sigma k\right)}}^{\Phi}|\textrm{ (tau-k)}$"];
                objP.Grid.Display         = [1];
                objP.Legend.Display       = [1];
                objP.Legend.FontSize      = [30.00];
                objP.Legend.Label         = ["$|\bar{\tau}_{f}^{\Phi}\hspace{0.800em}|$",...
                                             "$|\bar{\tau}_{f_{\left(\Sigma k\right)}}^{\Phi}|$",...
                                             "$|\bar{\tau}_{c}^{\Phi}|$",...
                                             "$|e_{c}^{\hspace{0.085em}\Phi}|$"];
                objP.Legend.Location      = ['NorthEast'];
                objP.Legend.NumColumns    = [2];
                objP.pbaspect             = [1,1,1];
                Zi                        = [cell(1,2)];
                %  > ------------------------------------------------------
                for i = [1,2]
                    %  > --------------------------------------------------
                    for j = 1:numel(Z)
                        Zi{i}(j,:) = Z{j}(i,:);
                    end
                    %  > --------------------------------------------------
                    objP.Axes.Label = ["",join(["$\left\|e\right\|_{_{",objP.S(i),"}}/\left\|e_{_{0}}\right\|_{_{",objP.S(i),"}}$"],'')];
                    %  > --------------------------------------------------
                    fprintf('Plotting...\n');
                    figure ('Color',objP.C(11,:),'Windowstate','Maximized');
                    hold on;
                    PP_1D.P3(objP,[1,5,2,3],Zi{i});
                    PP_1D.P4(objP);
                    if objP.E
                        exportgraphics(gcf,LF{3}(i),'ContentType','Vector');
                        movefile      (    LF{3}(i),PF);
                    end
                    %  > --------------------------------------------------
                end
                %  > ------------------------------------------------------
            end
            %  > ----------------------------------------------------------
            %  > #4.
            if any(F{4},2)
                %  > ------------------------------------------------------
                objP.Axes.Clipping        = [1];
                objP.Axes.FontSize{1}     = [25.00,25.00];
                objP.Axes.FontSize{2}     = [37.50,32.50];
                objP.Axes.Label           = ["$x$","$\textrm{Error magnitude}$"];
                objP.Axes.Reverse         = [0,0];
                objP.Axes.Scale           = [0,1];
                objP.Axes.Tick.X          = [0];
                objP.Axes.Tick.XTick      = [linspace(mX,MX,5)];
                objP.Axes.Tick.XTickLabel = [];
                objP.Axes.Tick.Y          = [0];
                objP.Axes.Tick.YTick      = [];
                objP.Axes.Tick.YTickLabel = [];
                objP.ET                   = [1e-10];
                objP.Grid.Display         = [1];
                objP.Legend.Display       = [1];
                objP.Legend.FontSize      = [30.00];
                objP.Legend.Location      = ['NorthEast'];
                objP.Legend.NumColumns    = [2];
                objP.LW                   = [3.50,2.50];
                objP.pbaspect             = [1,1,1];
                L                         = [cell(1,numel(F{4}))];
                M                         = [cell(1,numel(F{4}))];
                N                         = [cell(1,numel(F{4}))];
                X                         = [cell(1,numel(F{4}))];
                Y                         = [cell(1,numel(F{4}))];
                %  > ------------------------------------------------------
                for i = 1:numel(F{4})
                    if F{4}(i)
                        %  > ----------------------------------------------
                        switch i
                            case 1
                                L{i} = ["$|\bar{\tau}_{f_{\left(\Sigma k\right)}}^{\Phi}|$",...
                                        "$|e_{c}^{\hspace{0.085em}\Phi}\hspace{0.800em}|$",...
                                        "$|\bar{\tau}_{f_{\left(k\right)}}^{\Phi-\phi_{\left(h\right)}}|$",...
                                        "$|e_{c}^{\hspace{0.085em}\Phi-\phi_{\left(h\right)}}|$"];
                                M{i} = ["-o","-o";":","--"];
                                N{i} = [5,3;10,10];
                                X{i} = {msh.f.Xv,msh.c.Xc;
                                        msh.f.Xv,msh.c.Xc};
                                Y{i} = {obj(end).e.a   .t.f_abs(:,3),obj(end).e.a   .c.c_abs(:,1);
                                        obj(end).e.d(2).t.f_abs(:,3),obj(end).e.d(2).c.c_abs(:,1)};
                            case 2
                                L{i} = ["$|\bar{\tau}_{f_{\left(\Sigma k\right)}}^{\Phi}|$",...
                                        "$|e_{c}^{\hspace{0.085em}\Phi}\hspace{0.800em}|$",...
                                        "$|\bar{\tau}_{f_{\left(k\right)}}^{\Phi-\phi_{\left(l\right)}}|$",...
                                        "$|e_{c}^{\hspace{0.085em}\Phi-\phi_{\left(l\right)}}|$"];
                                M{i} = ["-o","-o";":","--"];
                                N{i} = [5,3;10,10];
                                X{i} = {msh.f.Xv,msh.c.Xc;
                                        msh.f.Xv,msh.c.Xc};
                                Y{i} = {obj(end).e.a   .t.f_abs(:,3),obj(end).e.a   .c.c_abs(:,1);
                                        obj(end).e.d(1).t.f_abs(:,3),obj(end).e.d(1).c.c_abs(:,1)};
                        end
                        %  > ----------------------------------------------
                        objP.Axes.Limits  = [PP_1D.AF_2(objP,X{i},[],Y{i},[])];
                        objP.Legend.Label = [L{i}];
                        %  > ----------------------------------------------
                        fprintf('Plotting...\n');
                        figure ('Color',objP.C(11,:),'Windowstate','Maximized');
                        hold on;
                        PP_1D.P1(...
                            [objP],...
                            [M{i}],...
                            [N{i}],...
                            [X{i}],...
                            [Y{i}]);
                        if objP.E
                            exportgraphics(gcf,LF{4}(i),'ContentType','Vector');
                            movefile      (    LF{4}(i),PF);
                        end
                        %  > ----------------------------------------------
                    end
                end
                %  > ------------------------------------------------------
            end
            %  > ----------------------------------------------------------
            %  > #5.
            if any(F{5},2)
                %  > ------------------------------------------------------
                objP.Axes.Clipping        = [1];
                objP.Axes.FontSize{1}     = [25.00,25.00];
                objP.Axes.FontSize{2}     = [37.50,32.50];
                objP.Axes.Label           = ["$x$","$\textrm{Error magnitude}$"];
                objP.Axes.Reverse         = [0,0];
                objP.Axes.Scale           = [0,1];
                objP.Axes.Tick.X          = [0];
                objP.Axes.Tick.XTick      = [linspace(mX,MX,5)];
                objP.Axes.Tick.XTickLabel = [];
                objP.Axes.Tick.Y          = [0];
                objP.Axes.Tick.YTick      = [];
                objP.Axes.Tick.YTickLabel = [];
                objP.ET                   = [1e-10];
                objP.Grid.Display         = [1];
                objP.Legend.Display       = [1];
                objP.Legend.FontSize      = [30.00];
                objP.Legend.Location      = ['NorthEast'];
                objP.Legend.NumColumns    = [2];
                objP.LW                   = [3.50,2.50];
                objP.pbaspect             = [1,1,1];
                L                         = [cell(1,numel(F{5}))];
                M                         = [cell(1,numel(F{5}))];
                N                         = [cell(1,numel(F{5}))];
                X                         = [cell(1,numel(F{5}))];
                Y                         = [cell(1,numel(F{5}))];
                %  > ------------------------------------------------------
                for i = 1:numel(F{5})
                    if F{5}(i)
                        %  > ----------------------------------------------
                        switch i
                            case 1
                                L{i} = ["$|\bar{\tau}_{f_{\left(\Sigma k\right)}}^{\Phi}|$",...
                                        "$|e_{c}^{\hspace{0.085em}\Phi}\hspace{0.800em}|$",...
                                        "$|\bar{\tau}_{f_{\left(k\right)}}^{\phi_{\left(h\right)}}|$",...
                                        "$|e_{c}^{\hspace{0.085em}\phi_{\left(h\right)}}|$"];
                                M{i} = ["-o","-o";":","--"];
                                N{i} = [5,3;10,10];
                                X{i} = {msh.f.Xv,msh.c.Xc;
                                        msh.f.Xv,msh.c.Xc};
                                Y{i} = {obj(end).e.a   .t.f_abs(:,3),obj(end).e.a   .c.c_abs(:,1);
                                        obj(end).e.x(2).t.f_abs(:,3),obj(end).e.x(2).c.c_abs(:,1)};
                            case 2
                                L{i} = ["$|\bar{\tau}_{f_{\left(\Sigma k\right)}}^{\Phi}|$",...
                                        "$|e_{c}^{\hspace{0.085em}\Phi}\hspace{0.800em}|$",...
                                        "$|\bar{\tau}_{f_{\left(k\right)}}^{\phi_{\left(l\right)}}|$",...
                                        "$|e_{c}^{\hspace{0.085em}\phi_{\left(l\right)}}|$"];
                                M{i} = ["-o","-o";":","--"];
                                N{i} = [5,3;10,10];
                                X{i} = {msh.f.Xv,msh.c.Xc;
                                        msh.f.Xv,msh.c.Xc};
                                Y{i} = {obj(end).e.a   .t.f_abs(:,3),obj(end).e.a   .c.c_abs(:,1);
                                        obj(end).e.x(1).t.f_abs(:,1),obj(end).e.x(1).c.c_abs(:,1)};
                        end
                        %  > ----------------------------------------------
                        objP.Axes.Limits  = [PP_1D.AF_2(objP,X{i},[],Y{i},[])];
                        objP.Legend.Label = [L{i}];
                        %  > ----------------------------------------------
                        fprintf('Plotting...\n');
                        figure ('Color',objP.C(11,:),'Windowstate','Maximized');
                        hold on;
                        PP_1D.P1(...
                            [objP],...
                            [M{i}],...
                            [N{i}],...
                            [X{i}],...
                            [Y{i}]);
                        if objP.E
                            exportgraphics(gcf,LF{5}(i),'ContentType','Vector');
                            movefile      (    LF{5}(i),PF);
                        end
                        %  > ----------------------------------------------
                    end
                end
            end
            %  > ----------------------------------------------------------
            %  > #6.
            if F{6}
                %  > ------------------------------------------------------
                switch inp.T
                    case 1
                        rl = unique(ceil(obj(end).s.u.p.c{1}./2)',"rows")';
                    case 2
                        for i = 1:numel(obj)
                            rl(:,i) = unique(ceil(obj{i}(end).s.u.p.c{1}./2)',"rows")'; %#ok<AGROW>
                        end
                    otherwise
                        return;
                end
                LA = ["$\Phi$",...
                      "$\phi_{\left(h\right)}$",...
                      "$\phi_{\left(\hspace{0.085em}l\hspace{0.085em}\right)}$"];
                MA = ["-","-",":"];
                NA = [1,2,3];
                %  > ------------------------------------------------------
                objP.Axes.Clipping         = [1];
                objP.Axes.FontSize{1}      = [25.00,25.00].*1.75;
                objP.Axes.FontSize{2}      = [37.50,32.50].*1.75;
                objP.Axes.Label            = ["$x$","$\textrm{Refinement level}$"];
                objP.Axes.Limits.X         = [min(msh.f.Xv(:,1),[],1),max(msh.f.Xv(:,1),[],1)];
                if max(rl,[],'all') == 1
                    objP.Axes.Limits.Y     = [1,max(rl,[],'all')+1];
                else
                    objP.Axes.Limits.Y     = [1,max(rl,[],'all')];
                end
                objP.Axes.Reverse          = [0,0];
                objP.Axes.Scale            = [0,0];
                objP.Axes.Tick.X           = [0];
                objP.Axes.Tick.XTick       = [linspace(mX,MX,5)];
                objP.Axes.Tick.XTickLabel  = [];
                objP.Axes.Tick.Y           = [0];
                objP.Axes.Tick.YTick       = [objP.Axes.Limits.Y(1):objP.Axes.Limits.Y(2)];
                objP.Axes.Tick.YTickLabel  = [];
                objP.Grid.Display          = [1];
                if inp.T == 1 || (inp.T == 2 && numel(obj) == 1)
                    objP.Legend.Display    = [0];
                    MB                     = [MA];
                    NB                     = [NA];
                else
                    MB                     = [MA(Z{:})];
                    NB                     = [NA(Z{:})];
                    objP.Legend.Display    = [1];
                    objP.Legend.FontSize   = [30.00].*1.75;
                    objP.Legend.Label      = [LA(Z{:})];
                    objP.Legend.Location   = ['NorthEast'];
                    objP.Legend.NumColumns = [1];
                end
                objP.LW                    = [10.00,5.00,5.00];
                objP.pbaspect              = [1,0.5,1];
                %  > ------------------------------------------------------
                fprintf('Plotting...\n');
                figure ('Color',objP.C(11,:),'Windowstate','Maximized');
                hold on;
                PP_1D.P5(...
                    [objP],...
                    [MB],...
                    [NB],...
                    [objP.Axes.Limits.X(1);msh.c.Xc;objP.Axes.Limits.X(2)],...
                    [rl;rl(end,:)]);
                if objP.E
                    exportgraphics(gcf,LF{6},'ContentType','Vector');
                    movefile      (    LF{6},PF);
                end
                %  > ------------------------------------------------------
            end
            %  > ----------------------------------------------------------
            %  > #7.
            if F{7}
                %  > ------------------------------------------------------
                h     = exp(1).^(linspace(log(1.00e-01),log(1.00e-03),10));
                n     = 1:10;
                dx    = cell(numel(n),numel(h));
                dX    = cell(numel(n),numel(h));
                cond2 = cell(2       ,numel(n));
                O     = cell(2       ,numel(n));
                for i = 1:numel(n)
                    for j = 1:numel(h)
                        dx   {i,j}      = ((-n(i)+1/2:1:n(i)-1/2).*h(j))';
                        dX   {i,j}      = (dx{i,j}-mean(dx{i,j},1))./std(dx{i,j},0,1);
                        cond2{1,i}(j,1) = cond(dx{i,j}.^(0:2*n(i)-1),2);
                        cond2{2,i}(j,1) = cond(dX{i,j}.^(0:2*n(i)-1),2);
                    end
                    O{1,i} = join(["$\mathbf{p}_{"      ,num2str(2.*i-1),"}$"],'');
                    O{2,i} = join(["$\mathbf{\hat{p}}_{",num2str(2.*i-1),"}$"],'');
                end
                %  > ------------------------------------------------------
                objP.Axes.Clipping        = [1];
                objP.Axes.FontSize{1}     = [25.00,25.00];
                objP.Axes.FontSize{2}     = [32.50,37.50];
                objP.Axes.Label           = ["$h$","$\kappa$"];
                objP.Axes.Limits.X        = [min(h,[],2),max(h,[],2)];
                objP.Axes.Limits.Y        = [1,10.^16];
                objP.Axes.Reverse         = [1,0];
                objP.Axes.Scale           = [1,1];
                objP.Axes.Tick.X          = [0];
                objP.Axes.Tick.XTick      = [];
                objP.Axes.Tick.XTickLabel = [];
                objP.Axes.Tick.Y          = [1];
                objP.Axes.Tick.YTick      = [];
                objP.Axes.Tick.YTickLabel = [];
                objP.Grid.Display         = [1];
                objP.Legend.Display       = [1];
                objP.Legend.FontSize      = [30.00];
                objP.Legend.Location      = ['NorthEast'];
                objP.Legend.NumColumns    = [2];
                objP.LW                   = [3.50];
                objP.pbaspect             = [1,1,1];
                %  > ------------------------------------------------------
                for i = 1:size(cond2,1)
                    %  > --------------------------------------------------
                    objP.Legend.Label = [O(i,:)];
                    %  > --------------------------------------------------
                    if ~objP.E
                        if i == 1
                            fprintf('Plotting...\n');
                            figure ('Color',objP.C(11,:),'Windowstate','Maximized');
                        end
                        subplot(1,2,i);
                        hold on;
                    else
                        fprintf('Plotting...\n');
                        figure ('Color',objP.C(11,:),'Windowstate','Maximized');
                        hold on;
                    end
                    plot(h,repmat(1./eps,1,numel(h)),"--",'Color',objP.C(10,:),'LineWidth',0.50,'MarkerFaceColor',objP.C(1,:));
                    PP_1D.P1(...
                        [objP],...
                        [repmat("-" ,1,numel(n))],...
                        [1:numel(n)],...
                        [repmat({h'},1,numel(n))],...
                        [cond2(i,:)]);
                    if objP.E
                        exportgraphics(gcf,LF{7}(i),'ContentType','Vector');
                        movefile      (    LF{7}(i),PF);
                    end
                    %  > --------------------------------------------------
                end
                %  > ------------------------------------------------------
            end
            %  > ----------------------------------------------------------
            %  > #8.
            if any(F{8},2)
                %  > ------------------------------------------------------
                h = [reshape(arrayfun(@(x) x.h,msh),[],1)];
                H = [h./h(1)];
                %  > ------------------------------------------------------
                objP.Axes.Clipping        = [1];
                objP.Axes.FontSize{1}     = [25.00,25.00];
                objP.Axes.FontSize{2}     = [32.50,32.50];
                objP.Axes.Label           = ["$h/h_{0}$","$\textrm{Error magnitude}$"];
                objP.Axes.Reverse         = [1,0];
                objP.Axes.Scale           = [1,1];
                objP.Axes.Tick.X          = [0];
                objP.Axes.Tick.XTick      = [H(end),H(1)];
                objP.Axes.Tick.XTickLabel = [];
                objP.Axes.Tick.Y          = [0];
                objP.Axes.Tick.YTick      = [];
                objP.Axes.Tick.YTickLabel = [];
                objP.ET                   = [1e-10];
                objP.Grid.Display         = [1];
                objP.Legend.Display       = [1];
                objP.Legend.FontSize      = [30.00];
                objP.Legend.Location      = ['NorthEast'];
                objP.Legend.NumColumns    = [2];
                objP.LW                   = [5.00,2.50];
                objP.pbaspect             = [1,1,1];
                L                         = [cell(1,numel(F{8}))];
                N                         = [cell(1,numel(F{8}))];
                P                         = [cell(numel(F{8}),2)];
                Y                         = [cell(numel(F{8}),2)];
                Z                         = [cell(numel(F{8}),2)];
                %  > ------------------------------------------------------
                switch ms
                    case 1
                        H0 = [1e-01];
                        L0 = [1e-09,1e+01];
                    case 2
                        H0 = [1e+02];
                        L0 = [1e-06,1e+04];
                end
                %  > ------------------------------------------------------
                for i = 1:numel(F{8})
                    if F{8}(i)
                        for j = [1,2]
                            %  > ------------------------------------------
                            switch i
                                case 1
                                    L{i}   = [join(["$\|\bar{\tau}_{f}^{\phi_{\left(h\right)}}\|_{_{",objP.S(j),"}}$"],''),...
                                              join(["$\|e_{c}^{\hspace{0.085em}\phi_{\left(h\right)}}\|_{_{",objP.S(j),"}}$"],''),...
                                              join(["$\|\bar{\tau}_{f}^{\Phi-\phi_{\left(h\right)}}\|_{_{",objP.S(j),"}}$"],''),...
                                              join(["$\|e_{c}^{\hspace{0.085em}\Phi-\phi_{\left(h\right)}}\|_{_{",objP.S(j),"}}$"],'')];
                                    N{i}   = [1,3;1,3];
                                    P{i,j} = ["2","4"];
                                    Y{i,j} = {arrayfun(@(x) x.e.x(2).n_abs.t.f(j,1),obj),arrayfun(@(x) x.e.x(2).n_abs.c(j,1),obj);...
                                              arrayfun(@(x) x.e.d(2).n_abs.t.f(j,1),obj),arrayfun(@(x) x.e.d(2).n_abs.c(j,1),obj)};
                                    Z{i,j} = {@(x) H0.*(1./H(1).*x).^2,@(x) H0.*(1./H(1).*x).^4};
                                case 2
                                    L{i}   = [join(["$\|\bar{\tau}_{f_{\left(0\right)}}^{\phi_{\left(h\right)}}\|_{_{",objP.S(j),"}}$"]),...
                                              join(["$\|\bar{\tau}_{f_{\left(1\right)}}^{\phi_{\left(h\right)}}\|_{_{",objP.S(j),"}}$"]),...
                                              join(["$\|\bar{\tau}_{f_{\left(0\right)}}^{\Phi-\phi_{\left(h\right)}}\|_{_{",objP.S(j),"}}$"]),...
                                              join(["$\|\bar{\tau}_{f_{\left(1\right)}}^{\Phi-\phi_{\left(h\right)}}+\ldots\|_{_{",objP.S(j),"}}$"])];
                                    N{i}   = [4,5;4,5];
                                    P{i,j} = ["2","4"];
                                    Y{i,j} = {arrayfun(@(x) x.e.x(2).n_abs.t.f(j,2),obj),arrayfun(@(x) x.e.x(2).n_abs.t.f(j,3),obj);...
                                              arrayfun(@(x) x.e.d(2).n_abs.t.f(j,2),obj),arrayfun(@(x) x.e.d(2).n_abs.t.f(j,3),obj)};
                                    Z{i,j} = {@(x) H0.*(1./H(1).*x).^2,@(x) H0.*(1./H(1).*x).^4};
                                case 3
                                    L{i}   = [join(["$\|\bar{\tau}_{f}^{\phi_{\left(l\right)}}\|_{_{",objP.S(j),"}}$"],''),...
                                              join(["$\|e_{c}^{\hspace{0.085em}\phi_{\left(l\right)}}\|_{_{",objP.S(j),"}}$"],''),...
                                              join(["$\|\bar{\tau}_{f}^{\Phi-\phi_{\left(l\right)}}\|_{_{",objP.S(j),"}}$"],''),...
                                              join(["$\|e_{c}^{\hspace{0.085em}\Phi-\phi_{\left(l\right)}}\|_{_{",objP.S(j),"}}$"],'')];
                                    N{i}   = [1,3;1,3];
                                    Y{i,j} = {arrayfun(@(x) x.e.x(1).n_abs.t.f(j,1),obj),arrayfun(@(x) x.e.x(1).n_abs.c(j,1),obj);...
                                              arrayfun(@(x) x.e.d(1).n_abs.t.f(j,1),obj),arrayfun(@(x) x.e.d(1).n_abs.c(j,1),obj)};
                                    switch j
                                        case 1
                                            P{i,j} = ["2","4"];
                                            Z{i,j} = {@(x) H0.*(1./H(1).*x).^2,@(x) H0.*(1./H(1).*x).^4};
                                        case 2
                                            switch ms
                                                case 1
                                                    P{i,j} = ["2","3"];
                                                    Z{i,j} = {@(x) H0.*(1./H(1).*x).^2,@(x) H0.*(1./H(1).*x).^3};
                                                case 2
                                                    P{i,j} = ["2","4"];
                                                    Z{i,j} = {@(x) H0.*(1./H(1).*x).^2,@(x) H0.*(1./H(1).*x).^4};
                                                otherwise
                                                    return;
                                            end
                                        otherwise
                                            return;
                                    end
                                case 4
                                    L{i}   = [join(["$\|\bar{\tau}_{f_{\left(0\right)}}^{\Phi}\hspace{0.095em}\|_{_{",objP.S(j),"}}$"]),...
                                              join(["$\|\bar{\tau}_{f_{\left(1\right)}}^{\phi_{\left(l\right)}}\|_{_{",objP.S(j),"}}$"]),...
                                              "",...
                                              join(["$\|\bar{\tau}_{f_{\left(1\right)}}^{\Phi-\phi_{\left(l\right)}}+\ldots\|_{_{",objP.S(j),"}}$"])];
                                    N{i}   = [4,5;11,5];
                                    Y{i,j} = {arrayfun(@(x) x.e.a.n_abs.t.f(j,2),obj),arrayfun(@(x) x.e.x(1).n_abs.t.f(j,1),obj);...
                                              NaN(size(H,1),size(H,2))               ,arrayfun(@(x) x.e.d(1).n_abs.t.f(j,3),obj)};
                                    switch j
                                        case 1
                                            switch ms
                                                case 1
                                                    P{i,j} = ["2","3"];
                                                    Z{i,j} = {@(x) H0.*(1./H(1).*x).^2,@(x) H0.*(1./H(1).*x).^3};
                                                case 2
                                                    P{i,j} = ["2","4"];
                                                    Z{i,j} = {@(x) H0.*(1./H(1).*x).^2,@(x) H0.*(1./H(1).*x).^4};
                                                otherwise
                                                    return;
                                            end
                                        case 2
                                            switch ms
                                                case 1
                                                    P{i,j} = ["2"];
                                                    Z{i,j} = {@(x) H0.*(1./H(1).*x).^2};
                                                case 2
                                                    P{i,j} = ["2","4"];
                                                    Z{i,j} = {@(x) H0.*(1./H(1).*x).^2,@(x) H0.*(1./H(1).*x).^4};
                                                otherwise
                                                    return;
                                            end
                                        otherwise
                                            return;
                                    end
                            end
                            %  > ------------------------------------------
                            objP.Axes.Limits  = [PP_1D.AF_2(objP,H,[],Y{i,j},[L0])];
                            objP.Legend.Label = [L{i}];
                            %  > ------------------------------------------
                            if false
                                S = AuxiliaryFunctions.Slope(H,cat(2,Y{i,j}{:})); disp(S); %#ok<UNRCH>
                            end
                            %  > ------------------------------------------
                            fprintf('Plotting...\n');
                            figure ('Color',objP.C(11,:),'Windowstate','Maximized');
                            hold on;
                            PP_1D.P1(objP,...
                                ["-o","-o";"-","-"],...
                                [N{i}],...
                                [repmat({H},size(Y{i,j}))],...
                                [Y{i,j}]);
                            PP_1D.P6(objP,...
                                [h],...
                                [8,12],...
                                [0.90,0.60],...
                                [P{i,j}],...
                                [Z{i,j}]);
                            if objP.E
                                exportgraphics(gcf,LF{8}(i,j),'ContentType','Vector');
                                movefile      (    LF{8}(i,j),PF);
                            end
                            %  > ------------------------------------------
                        end
                    end
                end
                %  > ------------------------------------------------------
            end
            %  > ----------------------------------------------------------
        end
        % >> 1.2. ---------------------------------------------------------
        %  > 1.2.1. -------------------------------------------------------
        %  > Auxiliary function #1 (fig. #1/#4/#5/#7/#8).
        function [] = P1(objP,M,N,X,Y)
            P = cell(size(Y));
            for i = 1:size(Y,1)
                for j = 1:size(Y,2)
                    P{i}{j} = plot(X{i,j},Y{i,j},M(i,j),'Color',objP.C(N(i,j),:),'LineWidth',objP.LW(i),'MarkerFaceColor',objP.C(N(i,j),:),'MarkerSize',objP.MS);
                end
            end
            PP_1D.AF_1(objP,cat(2,P{:}));
        end
        %  > 1.2.2. -------------------------------------------------------
        %  > Auxiliary function #2 (fig. #2).
        function [] = P2(objP,M,N,X,Y,Z)
            P = cell(size(Y));
            for i = 1:size(Z,1)
                for j = 1:size(Z,2)
                    if ~isempty(Z{i,j})
                        for m = 1:size(Z{i,j},1)
                            patch([Z{i,j}(m,1),Z{i,j}(m,2),Z{i,j}(m,2),Z{i,j}(m,1)],[objP.Axes.Limits.Y(1),objP.Axes.Limits.Y(1),objP.Axes.Limits.Y(2),objP.Axes.Limits.Y(2)],objP.C(N(1,j),:),'FaceAlpha',0.25,'Linestyle','None');
                        end
                    end
                end
            end
            for i = 1:size(Y,1)
                for j = 1:size(Y,2)
                    P{i}{j} = plot(X{i,j},Y{i,j},M(i,j),'Color',objP.C(N(2,j),:),'LineWidth',objP.LW(i),'MarkerFaceColor',objP.C(N(2,j),:),'MarkerSize',objP.MS);
                end
            end
            PP_1D.AF_1(objP,cat(2,P{:}));
        end
        %  > 1.2.3. -------------------------------------------------------
        %  > Auxiliary function #3 (fig. #3).
        function [] = P3(objP,N,Z)
            %  > Auxiliary variables.
            dX = [0.70];
            X1 = [objP.Axes.Limits.X(1)];
            X2 = [1+dX.*[-0.375,-0.125,0.125,0.375];
                  2+dX.*[-0.375,-0.125,0.125,0.375]];
            for i = 1:size(Z,1)
                for j = 1:size(Z,2)
                    plot([X1,X2(i,j)],[Z(i,j),Z(i,j)],'-.','Color',objP.C(N(j),:),'LineWidth',0.50,'MarkerFaceColor',objP.C(N(j),:));
                end
            end
            B = bar(Z,'BarWidth',0.85,'FaceColor','Flat');
            for i = 1:size(B,2)
                B(i).CData = objP.C(N(i),:);
            end
            PP_1D.AF_1(objP,{B});
            set(get(gca,'XAxis'),'TickLength',[0,0]);
            set(gca,'XGrid','Off');
        end
        %  > 1.2.4. -------------------------------------------------------
        %  > Auxiliary function #4 (fig. #3).
        function [] = P4(objP)
            %  > Auxiliary variables.
            d  = 0.1;
            Yb = objP.Axes.Tick.YTick(end)-d;
            
            if Yb <= 1-d
                patch([objP.Axes.Limits.X(1)-d./2,objP.Axes.Limits.X(1)+d./2,...
                       objP.Axes.Limits.X(1)+d./2,objP.Axes.Limits.X(1)-d./2],...
                      [Yb-d./4,Yb-d./4,...
                       Yb+d./4,Yb+d./4],objP.C(11,:),'EdgeColor',objP.C(11,:));
                plot ([objP.Axes.Limits.X(1)     ,objP.Axes.Limits.X(1)-d./2;
                       objP.Axes.Limits.X(1)-d./2,objP.Axes.Limits.X(1)+d./2;
                       objP.Axes.Limits.X(1)+d./2,objP.Axes.Limits.X(1)],...
                      [Yb-d./4,Yb-d./8;
                       Yb-d./8,Yb+d./8;
                       Yb+d./8,Yb+d./4],'k','LineWidth',0.5);
            end
        end
        %  > 1.2.5. -------------------------------------------------------
        %  > Auxiliary function #5 (fig. #6).
        %  > Add 1.15 inches to the left margin...
        function [] = P5(objP,M,N,X,Y)
            P = cell(1,size(Y,2));
            for i = 1:size(Y,2)
                P{i} = stairs(X,Y(:,i),M(i),'Color',objP.C(N(i),:),'Linewidth',objP.LW(i));
            end
            PP_1D.AF_1(objP,P);
        end
        %  > 1.2.6. -------------------------------------------------------
        %  > Auxiliary function #6 (fig. #7).
        function [] = P6(objP,h,x,y,P,Z)
            %  > Auxiliary variables.
            A = [exp(1).^(linspace(log(1),log(h(end)./h(1)),20))];
            B = [A(x(2)),A(x(1))];
            
            for i = 1:numel(Z)
                fplot(@(x) Z{i}(x),'-','Color',objP.C(10,:),'LineWidth',0.50);
                 plot([B(1),B(1);
                       B(1),B(1);
                       B(2),B(1)],...
                      [Z{i}(B(1)),Z{i}(B(2));
                       Z{i}(B(2)),Z{i}(B(2));
                       Z{i}(B(2)),Z{i}(B(1))],'Color',objP.C(10,:),'LineWidth',0.50);
                 text(min(B)-5e-03,mean(Z{i}(B).*y(i),2),P(i),'FontSize',22.50);
            end
        end
        % >> 1.3. ---------------------------------------------------------
        %  > 1.3.1. -------------------------------------------------------
        %  > Auxiliary function #1.
        function [] = AF_1(objP,P)
            %  > Axes.
            str.X = ["XTick","XTickLabel"];
            str.Y = ["YTick","YTickLabel"];
            for i = ["X","Y"]
                if isempty(objP.Axes.Tick.(str.(i)(2)))
                    if ~isempty(objP.Axes.Tick.(str.(i)(1)))
                        objP.Axes.Tick.(str.(i)(2)) = cell(numel(objP.Axes.Tick.(str.(i)(1))),1);
                        for j = 1:size(objP.Axes.Tick.(str.(i)(2)),1)
                            objP.Axes.Tick.(str.(i)(2)){j} = join(["$",num2str(objP.Axes.Tick.(str.(i)(1))(j)),"$"],'');
                        end
                    else
                        Lim                         = round(log10(objP.Axes.Limits.(i)),1);
                        objP.Axes.Tick.(str.(i)(2)) = cell (diff(Lim,1)+1,1);
                        for j = 1:size(objP.Axes.Tick.(str.(i)(2)),1)
                            if Lim(1)-1+j < 0
                                if numel(num2str(abs(Lim(1)-1+j))) == 1
                                    objP.Axes.Tick.(str.(i)(2)){j} = join(["$10^{-0",num2str(abs(Lim(1)-1+j)),"}$"],'');
                                else
                                    objP.Axes.Tick.(str.(i)(2)){j} = join(["$10^{- ",num2str(abs(Lim(1)-1+j)),"}$"],'');
                                end
                            else
                                if numel(num2str(abs(Lim(1)-1+j))) == 1
                                    objP.Axes.Tick.(str.(i)(2)){j} = join(["$10^{+0",num2str(abs(Lim(1)-1+j)),"}$"],'');
                                else
                                    objP.Axes.Tick.(str.(i)(2)){j} = join(["$10^{+ ",num2str(abs(Lim(1)-1+j)),"}$"],'');
                                end
                            end
                            objP.Axes.Tick.(str.(i)(1))(j) = 10.^(Lim(1)-1+j);
                        end
                        if ~objP.Axes.Tick.(i)
                            nrmv.(i) = true(size (objP.Axes.Tick.(str.(i)(1))));
                        else
                            nrmv.(i) = ~mod(log10(objP.Axes.Tick.(str.(i)(1))),2);
                        end
                        for j = 1:numel(str.(i))
                            objP.Axes.Tick.(str.(i)(j)) = objP.Axes.Tick.(str.(i)(j))(nrmv.(i));
                        end
                    end
                end
            end
            set(gca,'Box','On');
            if objP.Axes.Clipping
                set(gca,'Clipping','On');
            end
            if objP.Axes.Reverse(1)
                set(gca,'XDir','Reverse');
            end
            if objP.Axes.Reverse(2)
                set(gca,'YDir','Reverse');
            end
            if objP.Axes.Scale(1)
                set(gca,'XScale','Log');
            end
            if objP.Axes.Scale(2)
                set(gca,'YScale','Log');
            end
            set(gca,'XLim',objP.Axes.Limits.X,'XTick',objP.Axes.Tick.XTick,'XTickLabel',objP.Axes.Tick.XTickLabel); set(get(gca,'XAxis'),'Fontsize',objP.Axes.FontSize{1}(1)); xlabel(objP.Axes.Label(1),'Fontsize',objP.Axes.FontSize{2}(1)); xtickangle(0);
            set(gca,'YLim',objP.Axes.Limits.Y,'YTick',objP.Axes.Tick.YTick,'YTickLabel',objP.Axes.Tick.YTickLabel); set(get(gca,'YAxis'),'Fontsize',objP.Axes.FontSize{1}(2)); ylabel(objP.Axes.Label(2),'Fontsize',objP.Axes.FontSize{2}(2)); ytickangle(0);
            pbaspect(objP.pbaspect);
            %  > Grid.
            if objP.Grid.Display
                grid(gca,'On'); set(gca,'GridLineStyle','-','GridAlpha',0.15,'MinorGridAlpha',0.00);
            end
            %  > Interpreter.
            fn = fieldnames(get(groot,'factory'));
            fc = find(contains(fn,'Interpreter'));
            for i = 1:numel(fc)
                set(groot,strrep(fn{fc(i)},'factory','Default'),'Latex');
            end
            %  > Legend.
            if objP.Legend.Display
                legend([P{:}],objP.Legend.Label,'AutoUpdate','Off','Location',objP.Legend.Location,'FontSize',objP.Legend.FontSize,'NumColumns',objP.Legend.NumColumns);
            end
        end
        %  > 1.3.2. -------------------------------------------------------
        %  > Auxiliary function #2.
        function [Limits] = AF_2(objP,V1,V2,V3,V4)
            str.X = ["XTick"];
            str.Y = ["YTick"];
            for i = ["X","Y"]
                switch i
                    case "X", V.(i) = {V1,V2}; 
                    case "Y", V.(i) = {V3,V4};
                end
                if ~isempty(V.(i){2})
                    Limits.(i) = V.(i){2};
                else
                    if ~isempty(objP.Axes.Tick.(str.(i)))
                        Limits.(i) = [min(objP.Axes.Tick.(str.(i)),[],2),max(objP.Axes.Tick.(str.(i)),[],2)];
                    else
                        Vm = reshape(cellfun(@min,V.(i){1}),1,[]); Vm = Vm(~isinf(Vm) & ~isnan(Vm)); Limits.(i)(1) = min(Vm,[],2);
                        VM = reshape(cellfun(@max,V.(i){1}),1,[]); VM = VM(~isinf(VM) & ~isnan(VM)); Limits.(i)(2) = max(VM,[],2);
                        if  Limits.(i)(1) < objP.ET
                            Limits.(i)(1) = objP.ET;
                            Limits.(i)(2) = 10.^(ceil(log10(Limits.(i)(2)))+0);
                        else
                            Limits.(i)(1) = 10.^(ceil(log10(Limits.(i)(1)))-1);
                            Limits.(i)(2) = 10.^(ceil(log10(Limits.(i)(2)))+0);
                        end
                        if Limits.(i)(1) == Limits.(i)(2)
                            Limits.(i)(1) = Limits.(i)(2)./10;
                        end
                    end
                end
            end
        end
    end
end