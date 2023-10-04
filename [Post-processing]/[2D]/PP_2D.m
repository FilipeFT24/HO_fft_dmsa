classdef PP_2D
    methods (Static)
        %% > 1. -----------------------------------------------------------
        % >> 1.1. ---------------------------------------------------------
        %  > 2D plot(s).
        function [] = Plot(inp,msh,obj,F)
            % >> ----------------------------------------------------------
            % >> Auxiliary variables.
            %  > ----------------------------------------------------------
            objP.C       = [linspecer(9,'Qualitative')]; %  > Color.
            objP.C(10,:) = [0,0,0];                      %  > k.
            objP.C(11,:) = [1,1,1];                      %  > w.
            objP.E       = [0];                          %  > Export?
            objP.S       = ["1","\infty"];               %  > Lk error norm(s).
            if size(msh,2) == 1
                X = cat(1,msh.f.xy.v{:}); mX = min(X,[],1);
                                          MX = max(X,[],1);
            end
            if ~inp.PS.EE
                F{3}(1:2) = 0;
            end
            if size(obj,1) == 1
                F{6}(1:3) = 0;
            else
                F{1}(1)   = 0;
                F{2}(1:2) = 0;
                F{3}(1:2) = 0;
                F{4}(1)   = 0;
                F{5}(1:2) = 0;
            end
            %  > ----------------------------------------------------------
            if objP.E
                LF = cell(numel(F),1);
                for i = 1:numel(F)
                    switch i
                        case 1
                            LF{i} = "DS.pdf";
                        case 2
                            for j = 1:numel(F{i})
                                LF{i}(j) = join(["NS (1.1.",num2str(j),").pdf"],'');
                            end
                        case 3
                            for j = 1:numel(F{i})
                                LF{i}(j) = join(["NS (1.2.",num2str(j),").pdf"],'');
                            end
                        case 4
                            LF{i} = "AMR (1.3).pdf";
                        case 5
                            for j = 1:numel(F{i})
                                for k = [1,2]
                                    LF{i}(j,k) = join(["B (",num2str(j),".",num2str(k),").pdf"],'');
                                end
                            end
                        case 6
                            A = [1,1;1,1;2,2];
                            B = [1,2;3,4;1,2];
                            for j = 1:numel(F{i})
                                for k = [1,2]
                                    LF{i}(j,k) = join(["NEE (D",num2str(A(j,k)),".",num2str(B(j,k)),").pdf"],'');
                                end
                            end
                    end
                end
            end
            PF = "[Post-processing]/[.PDF Files]/[2D]";
            %  > ----------------------------------------------------------
            
            % >> ----------------------------------------------------------
            % >> Plot...
            %  > ----------------------------------------------------------
            %  > #1.
            if any(F{1},2)
                %  > ------------------------------------------------------
                objP.Axes.FontSize{1} = [25.00,25.00];
                objP.Axes.FontSize{2} = [37.50,37.50];
                objP.Axes.Label       = ["$x$","$y$"];
                objP.ET               = [1e-10];
                objP.pbaspect         = [1,1,1];
                %  > ------------------------------------------------------
                objP.Axes.Limits.X        = [mX(1),MX(1)];
                objP.Axes.Limits.Y        = [mX(2),MX(2)];
                objP.Axes.Limits.Z        = [];
                objP.Axes.Tick.X          = [0];
                objP.Axes.Tick.XTick      = [linspace(mX(1),MX(1),5)];
                objP.Axes.Tick.XTickLabel = [];
                objP.Axes.Tick.Y          = [0];
                objP.Axes.Tick.YTick      = [linspace(mX(2),MX(2),5)];
                objP.Axes.Tick.YTickLabel = [];
                objP.Axes.Tick.Z          = [0];
                objP.Axes.Tick.ZTick      = [];
                objP.Axes.Tick.ZTickLabel = [];
                objP.Colorbar.Display     = [0];
                objP.LW                   = [5.00];
                objP.MS                   = [256.725.*msh.h.*(MX(1)-mX(1))./4];
                f                         = [1];
                %  > ------------------------------------------------------
                fprintf('Plotting...\n');
                figure ('Color',objP.C(11,:),'Windowstate','Maximized');
                hold on;
                cellfun(@(x) plot(x(:,1),x(:,2),'-','Color',objP.C(10,:),'Linewidth',0.50),msh.f.xy.v);
                PP_2D.P1(...
                    [objP],...
                    [msh],...
                    [f(1)],...
                    [obj(end).s.c.s.st{f(1),1}]);
                if objP.E
                    exportgraphics(gcf,LF{1}(i),'ContentType','Vector');
                    movefile      (    LF{1}(i),PF);
                end
                %  > ------------------------------------------------------
            end
            %  > ----------------------------------------------------------
            %  > #2.
            if any(F{2},2)
                %  > ------------------------------------------------------
                objP.Axes.FontSize{1}     = [25.00,25.00];
                objP.Axes.FontSize{2}     = [37.50,37.50];
                objP.Axes.Label           = ["$x$","$y$"];
                objP.Axes.Limits.X        = [mX(1),MX(1)];
                objP.Axes.Limits.Y        = [mX(2),MX(2)];
                objP.Axes.Limits.Z        = [];
                objP.Axes.Tick.X          = [0];
                objP.Axes.Tick.XTick      = [linspace(mX(1),MX(1),5)];
                objP.Axes.Tick.XTickLabel = [];
                objP.Axes.Tick.Y          = [0];
                objP.Axes.Tick.YTick      = [linspace(mX(2),MX(2),5)];
                objP.Axes.Tick.YTickLabel = [];
                objP.Axes.Tick.Z          = [0];
                objP.Axes.Tick.ZTick      = [];
                objP.Axes.Tick.ZTickLabel = [];
                objP.Colorbar.Display     = [1];
                objP.Colorbar.FontSize    = [25.00];
                objP.Colorbar.Map         = [1];
                objP.Colorbar.Position    = [-0.0750,-0.1150,0.0250,0.1150];
                objP.Colorbar.Scale       = [1];
                objP.ET                   = [1e-10];
                objP.pbaspect             = [1,1,1];
                Y                         = [cell(1,numel(F{2}))];
                %  > ------------------------------------------------------
                for i = 1:numel(F{2})
                    if F{2}(i)
                        %  > ----------------------------------------------
                        switch i
                            case 1, Y{i} = obj(end).e.a.t.c_abs(:,1);
                            case 2, Y{i} = obj(end).e.a.c.c_abs(:,1);
                        end
                        %  > ----------------------------------------------
                        fprintf('Plotting...\n');
                        figure ('Color',objP.C(11,:),'Windowstate','Maximized');
                        hold on;
                        PP_2D.P3(...
                            [objP],...
                            [msh.c.c.xy.v],...
                            [Y{i}]);
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
            if any(F{3},2)
                %  > ------------------------------------------------------
                dZ                        = {[1e-10,1e+00],...
                                             [1e-10,1e+00]};
                %  > ------------------------------------------------------
                objP.Axes.Clipping        = [1];
                objP.Axes.FontSize{1}     = [25.00,25.00];
                objP.Axes.FontSize{2}     = [32.50,32.50];
                objP.Axes.Reverse         = [0,0];
                objP.Axes.Scale           = [1,1];
                objP.Axes.Tick.X          = [1];
                objP.Axes.Tick.XTick      = [];
                objP.Axes.Tick.XTickLabel = [];
                objP.Axes.Tick.Y          = [1];
                objP.Axes.Tick.YTick      = [];
                objP.Axes.Tick.YTickLabel = [];
                objP.ET                   = [1e-10];
                objP.Grid.Display         = [1];
                objP.Legend.Display       = [1];
                objP.Legend.FontSize      = [30.00];
                objP.Legend.Location      = ['NorthWest'];
                objP.Legend.NumColumns    = [1];
                objP.LW                   = [3.50];
                objP.MS                   = [6.50];
                objP.pbaspect             = [1,1,1];
                Y                         = [cell(1,numel(F{3}))];
                %  > ------------------------------------------------------
                for i = 1:numel(F{3})
                    if F{3}(i)
                        %  > ----------------------------------------------
                        switch i
                            case 1
                                objP.Axes.Label   = ["$|\bar{\tau}_{f_{\left(\Sigma k\right)}}^{\Phi}|$",...
                                                     "$|\bar{\tau}_{f_{\left(k\right)}}^{\phi}|$"];
                                objP.Legend.Label = ["$|\bar{\tau}_{f_{\left(k\right)}}^{\phi_{\left(h\right)}}|$",...
                                                     "$|\bar{\tau}_{f_{\left(k\right)}}^{\phi_{\left(l\right)}}\hspace{0.125em}|$"];
                                Y{i}              = {obj(end).e.a.t.f_abs(:,3),obj(end).e.x(2).t.f_abs(:,3);
                                                     obj(end).e.a.t.f_abs(:,3),obj(end).e.x(1).t.f_abs(:,1)};
                            case 2
                                objP.Axes.Label   = ["$|e_{c}^{\Phi}|$",...
                                                     "$|e_{c}^{\phi}|$"];
                                objP.Legend.Label = ["$|e_{c}^{\phi_{\left(h\right)}}|$",...
                                                     "$|e_{c}^{\phi_{\left(l\right)}}\hspace{0.125em}|$"];
                                Y{i}              = {obj(end).e.a.c.c_abs(:,1),obj(end).e.x(2).c.c_abs(:,1);
                                                     obj(end).e.a.c.c_abs(:,1),obj(end).e.x(1).c.c_abs(:,1)};
                        end
                        %  > ----------------------------------------------
                        objP.Axes.Limits = [PP_1D.AF_2(objP,[],dZ{i},[],dZ{i})];
                        %  > ----------------------------------------------
                        fprintf('Plotting...\n');
                        figure ('Color',objP.C(11,:),'Windowstate','Maximized');
                        hold on;
                        plot(dZ{i},dZ{i},'-','Color',objP.C(10,:));
                        PP_1D.P1(...
                            [objP],...
                            ["o","o"],...
                            [1,2],...
                            [Y{i}(:,1)'],...
                            [Y{i}(:,2)']);
                        if objP.E
                            exportgraphics(gcf,LF{3}(i),'ContentType','Vector');
                            movefile      (    LF{3}(i),PF);
                        end
                        %  > ----------------------------------------------
                    end
                end
                %  > ------------------------------------------------------
            end
            %  > ----------------------------------------------------------
            %  > #4.
            if F{4}
                %  > ------------------------------------------------------
                rl = unique(ceil(obj(end).s.u.p.c{1}./2)','rows')';
                l  = false (size(rl));
                Pu = unique(rl,'rows');
                Nu = numel (Pu);
                %  > ------------------------------------------------------
                objP.Axes.FontSize{1}     = [25.00,25.00];
                objP.Axes.FontSize{2}     = [37.50,37.50];
                objP.Axes.Label           = ["$x$","$y$"];
                objP.Axes.Limits.X        = [mX(1),MX(1)];
                objP.Axes.Limits.Y        = [mX(2),MX(2)];
                objP.Axes.Limits.Z        = [0.5,0.5+Nu];
                objP.Axes.Tick.X          = [0];
                objP.Axes.Tick.XTick      = [linspace(mX(1),MX(1),5)];
                objP.Axes.Tick.XTickLabel = [];
                objP.Axes.Tick.Y          = [0];
                objP.Axes.Tick.YTick      = [linspace(mX(2),MX(2),5)];
                objP.Axes.Tick.YTickLabel = [];
                objP.Axes.Tick.Z          = [0];
                objP.Axes.Tick.ZTick      = [1:Nu];
                objP.Axes.Tick.ZTickLabel = [string(Pu)];
                objP.Colorbar.Display     = [1];
                objP.Colorbar.FontSize    = [25.00];
                objP.Colorbar.Map         = [0];
                objP.Colorbar.Position    = [-0.0750,-0.1150,0.0250,0.1150];
                objP.Colorbar.Scale       = [0];
                objP.ET                   = [1e-10];
                objP.MS                   = [256.725.*msh.h.*(MX(1)-mX(1))./1.925];
                objP.pbaspect             = [1,1,1];
                %  > ------------------------------------------------------
                fprintf('Plotting...\n');
                figure ('Color',objP.C(11,:),'Windowstate','Maximized');
                hold on;
                cellfun(@(x) plot(x(:,1),x(:,2),'-','Color',objP.C(10,:),'Linewidth',0.50),msh.f.xy.v);
                PP_2D.P4(...
                    [objP],...
                    [msh],...
                    [l],...
                    [rl]);
                if objP.E
                    exportgraphics(gcf,LF{4},'ContentType','Vector');
                    movefile      (    LF{4},PF);
                end
                %  > ------------------------------------------------------
            end
            %  > ----------------------------------------------------------
            %  > #5.
            if any(F{5},2)
                %  > ------------------------------------------------------
                h   = exp(1).^(linspace(log(1.00e-01),log(1.00e-03),10));
                n   = 1:10;
                m   = 1:10;
                q   = {[ 3],...
                       [ 5, 3],...
                       [ 7, 5, 3],...
                       [ 9, 7, 7, 5],...
                       [11, 9, 9, 7, 5],...
                       [13,11,11, 9, 9, 5],...
                       [15,13,13,13,11, 9, 5],...
                       [17,15,15,15,13,11, 9, 5],...
                       [19,17,17,17,15,15,13,11, 7],...
                       [21,19,19,19,17,17,15,13,11, 7]};
                c   = cell(1,n(end));
                di  = cell(1,n(end));
                dij = cell(  n(end));
                d   = cell(n(end),numel(h));
                dx  = cell(n(end),numel(h));
                dX  = cell(n(end),numel(h));
                for i = n
                    c{i} = B1_2D.PR1(repmat(2.*i-1,1,2));
                end
                for i = n
                    for j = 1:i
                        dij{j,i} = @(x) [[repmat(x*(-j+1/2),q{i}(j),1)],x*[-(q{i}(j)-1)/2:(q{i}(j)-1)/2]';
                                         [repmat(x*( j-1/2),q{i}(j),1)],x*[-(q{i}(j)-1)/2:(q{i}(j)-1)/2]'];
                        %        = @(x) [[repmat(-(2.*(i-j)+1).*x./2,1,2.*j+1);-j*x:x:j*x]';...
                        %                [repmat( (2.*(i-j)+1).*x./2,1,2.*j+1);-j*x:x:j*x]'];
                    end
                    di{i} = @(x) cellfun(@(f) f(x),dij(:,i),'UniformOutput',false);
                    for j = 1:numel(h)
                        d{i,j} = di{i}(h(j));
                    end
                    for j = 1:numel(h)
                        dx{i,j} = cat(1,d{i,j}{:});
                        dX{i,j} = (dx{i,j}-mean(dx{i,j},1))./std(dx{i,j},0,1);
                    end
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
                objP.Axes.Tick.XTick      = [];
                objP.Axes.Tick.XTickLabel = [];
                objP.Axes.Tick.YTick      = [];
                objP.Axes.Tick.YTickLabel = [];
                objP.Grid.Display         = [1];
                objP.Legend.Display       = [1];
                objP.Legend.FontSize      = [30.00];
                objP.Legend.Location      = ['NorthEast'];
                objP.Legend.NumColumns    = [2];
                objP.LW                   = [3.50];
                objP.MS                   = [6.50];
                objP.pbaspect             = [1,1,1];
                cond2                     = [cell(1,numel(F{5}))];
                O                         = [cell(1,numel(F{5}))];
                %  > ------------------------------------------------------
                for i = 1:numel(F{5})
                    if F{5}(i)
                        %  > ----------------------------------------------
                        switch i
                            case 1
                                for j = n
                                    for k = 1:numel(h)
                                        cond2{i}{1,j}(k,1) = cond(c{j}(1).c(1,:).*dx{j,k}(:,1).^c{j}(1).e(1,:).*dx{j,k}(:,2).^c{j}(1).e(2,:));
                                        cond2{i}{2,j}(k,1) = cond(c{j}(1).c(1,:).*dX{j,k}(:,1).^c{j}(1).e(1,:).*dX{j,k}(:,2).^c{j}(1).e(2,:));
                                    end
                                    O{1}(:,j) = {join(["$\mathbf{p}_{",num2str(2.*j-1),"}$"],''),join(["$\mathbf{\hat{p}}_{",num2str(2.*j-1),"}$"],'')};
                                end
                                objP.Axes.Limits.Y = [1,10.^16];
                                objP.Axes.Tick.X   = [0];
                                objP.Axes.Tick.Y   = [1];
                            case 2
                                for j = 1:numel(m)
                                    for k = 1:numel(h)
                                        cond2{i}{1,j}(k,1) = cond(sqrt(sqrt(sum(dx{1,k}.^2,2)).^(-j)).*...
                                            c{1}(1).c(1,:).*dx{1,k}(:,1).^c{1}(1).e(1,:).*dx{1,k}(:,2).^c{1}(1).e(2,:));
                                        cond2{i}{2,j}(k,1) = cond(sqrt(sqrt(sum((dx{1,k}./std(dx{1,k},0,1)).^2,2)).^(-j)).*...
                                            c{1}(1).c(1,:).*dX{1,k}(:,1).^c{1}(1).e(1,:).*dX{1,k}(:,2).^c{1}(1).e(2,:));
                                    end
                                    O{2}(:,j) = {join(["$1/d_{i}^{\,",num2str(j),"}$"],''),join(["$1/\hat{d}_{i}^{\,",num2str(j),"}$"],'')};
                                end
                                objP.Axes.Limits.Y = [1,10.^05];
                                objP.Axes.Tick.X   = [0];
                                objP.Axes.Tick.Y   = [0];
                        end
                        %  > ----------------------------------------------
                        for j = 1:size(cond2{i},1)
                            %  > ------------------------------------------
                            objP.Legend.Label = [O{i}(j,:)];
                            %  > ------------------------------------------
                            fprintf('Plotting...\n');
                            figure ('Color',objP.C(11,:),'Windowstate','Maximized');
                            hold on;
                            plot(h,repmat(1./eps,1,numel(h)),"--",'Color',objP.C(10,:),'LineWidth',0.50,'MarkerFaceColor',objP.C(1,:));
                            PP_1D.P1(...
                                [objP],...
                                [repmat("-" ,1,numel(m))],...
                                [1:numel(m)],...
                                [repmat({h'},1,numel(m))],...
                                [cond2{i}(j,:)]);
                            if objP.E
                                exportgraphics(gcf,LF{5}(i,j),'ContentType','Vector');
                                movefile      (    LF{5}(i,j),PF);
                            end
                            %  > ------------------------------------------
                        end
                        %  > ----------------------------------------------
                    end
                end
                %  > ------------------------------------------------------
            end
            %  > ----------------------------------------------------------
            %  > #6.
            if any(F{6},2)
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
                objP.MS                   = [6.50];
                objP.pbaspect             = [1,1,1];
                L                         = [cell(1,numel(F{6}))];
                N                         = [cell(1,numel(F{6}))];
                P                         = [cell(numel(F{6}),2)];
                O                         = [cell(numel(F{6}),2)];
                Y                         = [cell(numel(F{6}),2)];
                Z                         = [cell(numel(F{6}),2)];
                %  > ------------------------------------------------------
                for i = 1:numel(F{6})
                    if F{6}(i)
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
                                    switch j
                                        case 1
                                            C      = [0.90,0.80,0.65,0.55];
                                            O{i,j} = [1e-07,1e+01];
                                            P{i,j} = ["2","3","4","5"];
                                            Z{i,j} = {@(x) 1e-01.*(1./H(1).*x).^2,@(x) 1e-02.*(1./H(1).*x).^3,@(x) 1e-01.*(1./H(1).*x).^4,@(x) 1e-02.*(1./H(1).*x).^5};
                                        case 2
                                            C      = [0.90,0.65];
                                            O{i,j} = [1e-06,1e+02];
                                            P{i,j} = ["2","4"];
                                            Z{i,j} = {@(x) 1e+00.*(1./H(1).*x).^2,@(x) 1e+00.*(1./H(1).*x).^4};
                                    end
                                case 2
                                    L{i}   = [join(["$\|\bar{\tau}_{f_{\left(0\right)}}^{\phi_{\left(h\right)}}\|_{_{",objP.S(j),"}}$"]),...
                                              join(["$\|\bar{\tau}_{f_{\left(1\right)}}^{\phi_{\left(h\right)}}\|_{_{",objP.S(j),"}}$"]),...
                                              join(["$\|\bar{\tau}_{f_{\left(0\right)}}^{\Phi-\phi_{\left(h\right)}}\|_{_{",objP.S(j),"}}$"]),...
                                              join(["$\|\bar{\tau}_{f_{\left(1\right)}}^{\Phi-\phi_{\left(h\right)}}+\ldots\|_{_{",objP.S(j),"}}$"])];
                                    N{i}   = [4,5;4,5];
                                    Y{i,j} = {arrayfun(@(x) x.e.x(2).n_abs.t.f(j,2),obj),arrayfun(@(x) x.e.x(2).n_abs.t.f(j,3),obj);...
                                              arrayfun(@(x) x.e.d(2).n_abs.t.f(j,2),obj),arrayfun(@(x) x.e.d(2).n_abs.t.f(j,3),obj)};
                                    switch j
                                        case 1
                                            C      = [0.80,0.55];
                                            O{i,j} = [1e-07,1e+01];
                                            P{i,j} = ["3","5"];
                                            Z{i,j} = {@(x) 1e-02.*(1./H(1).*x).^3,@(x) 1e-02.*(1./H(1).*x).^5};
                                        case 2
                                            C      = [0.90,0.65];
                                            O{i,j} = [1e-06,1e+02];
                                            P{i,j} = ["2","4"];
                                            Z{i,j} = {@(x) 1e+00.*(1./H(1).*x).^2,@(x) 1e+00.*(1./H(1).*x).^4};
                                    end
                                case 3
                                    L{i}   = [join(["$\|e_{c}^{\hspace{0.085em}\phi_{\left(l\right)}}\|_{_{",objP.S(j),"}}$"],''),...
                                              join(["$\|\bar{\tau}_{f_{\left(1\right)}}^{\phi_{\left(l\right)}}\|_{_{",objP.S(j),"}}$"]),...
                                              join(["$\|e_{c}^{\hspace{0.085em}\Phi-\phi_{\left(l\right)}}\|_{_{",objP.S(j),"}}$"],''),...
                                              join(["$\|\bar{\tau}_{f_{\left(1\right)}}^{\Phi-\phi_{\left(l\right)}}\|_{_{",objP.S(j),"}}$"])];
                                    N{i}   = [3,5;3,5];
                                    Y{i,j} = {arrayfun(@(x) x.e.x(1).n_abs.c(j,1),obj),arrayfun(@(x) x.e.x(1).n_abs.t.f(j,1),obj);...
                                              arrayfun(@(x) x.e.d(1).n_abs.c(j,1),obj),arrayfun(@(x) x.e.d(1).n_abs.t.f(j,3),obj)};
                                    switch j
                                        case 1
                                            C      = [0.90,0.80];
                                            O{i,j} = [1e-07,1e+01];
                                            P{i,j} = ["2","3"];
                                            Z{i,j} = {@(x) 1e-01.*(1./H(1).*x).^2,@(x) 1e-02.*(1./H(1).*x).^3};
                                        case 2
                                            C      = [0.90];
                                            O{i,j} = [1e-06,1e+02];
                                            P{i,j} = ["2"];
                                            Z{i,j} = {@(x) 1e-01.*(1./H(1).*x).^2};
                                    end
                            end
                            %  > ------------------------------------------
                            objP.Axes.Limits  = [PP_1D.AF_2(objP,H,[],Y{i,j},[O{i,j}])];
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
                                [6,14],...
                                [C],...
                                [P{i,j}],...
                                [Z{i,j}]);
                            if objP.E
                                exportgraphics(gcf,LF{6}(i,j),'ContentType','Vector');
                                movefile      (    LF{6}(i,j),PF);
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
        %  > Auxiliary function #1 (fig. #1).
        function [] = P1(objP,msh,f,s)
            %  > Cell(s).
            cellfun(@(x) patch(x(:,1),x(:,2),objP.C(1,:),'FaceAlpha',0.25,'Linestyle','None'),msh.c.c.xy.v(s.i(s.l)));
            %  > Cell centroid(s).
            plot(msh.c.c.xy.c(:,1),msh.c.c.xy.c(:,2),'o','Color',objP.C(10,:),'MarkerFaceColor',objP.C(10,:),'MarkerSize',objP.MS./4);
            %  > Discretisation stencil.
            plot(s.X(:,1),s.X(:,2),'s','Color',objP.C(1,:),'MarkerFaceColor',objP.C(1,:),'MarkerSize',objP.MS);
            %  > Face(s).
            plot(msh.f.xy.v{f}(:,1),msh.f.xy.v{f(1)}(:,2),'-','Color',objP.C(1,:),'LineWidth',objP.LW(1));
            PP_2D.AF_1(objP);
        end
        %  > 1.2.2. -------------------------------------------------------
        %  > Auxiliary function #2 (fig. #1).
        function [] = P2(objP,X,Y)
            cellfun(@(x,y) patch(x(:,1),x(:,2),y,'LineWidth',1.00),X,num2cell(Y));
            PP_2D.AF_2(objP,Y);
        end
        %  > 1.2.3. -------------------------------------------------------
        %  > Auxiliary function #3 (fig. #2).
        function [] = P3(objP,X,Y)
            cellfun(@(x,y) patch(x(:,1),x(:,2),y,'Linestyle','None'),X,num2cell(Y));
            PP_2D.AF_2(objP,Y);
        end
        %  > 1.2.4. -------------------------------------------------------
        %  > Auxiliary function #4 (fig. #4).
        function [] = P4(objP,msh,X,Y)
            cellfun(@(x,y) plot(x(:,1),x(:,2),'o','Color',objP.C(y,:),'MarkerFaceColor',objP.C(11,:),'MarkerSize',objP.MS),num2cell(msh.f.xy.c( X,:),2),num2cell(Y( X,1),2));
            cellfun(@(x,y) plot(x(:,1),x(:,2),'o','Color',objP.C(y,:),'MarkerFaceColor',objP.C(y ,:),'MarkerSize',objP.MS),num2cell(msh.f.xy.c(~X,:),2),num2cell(Y(~X,1),2));
            PP_2D.AF_2(objP,Y);
        end
        % >> 1.3. ---------------------------------------------------------
        %  > 1.3.1. -------------------------------------------------------
        %  > Auxiliary function #1.
        function [] = AF_1(objP)
            %  > Axes.
            str.X = ["XTick","XTickLabel"];
            str.Y = ["YTick","YTickLabel"];
            for i = ["X","Y"]
                objP.Axes.Tick.(str.(i)(2)) = cell(numel(objP.Axes.Tick.(str.(i)(1))),1);
                for j = 1:size(objP.Axes.Tick.(str.(i)(2)),1)
                    objP.Axes.Tick.(str.(i)(2)){j} = join(["$",num2str(objP.Axes.Tick.(str.(i)(1))(j)),"$"],'');
                end
            end
            set(gca,'Box','Off','Clipping','Off');
            set(gca,'XLim',objP.Axes.Limits.X,'XTick',objP.Axes.Tick.XTick,'XTickLabel',objP.Axes.Tick.XTickLabel,'TickLength',[0,0]); set(get(gca,'XAxis'),'Fontsize',objP.Axes.FontSize{1}(1)); xlabel(objP.Axes.Label(1),'Fontsize',objP.Axes.FontSize{2}(1)); xtickangle(0);
            set(gca,'YLim',objP.Axes.Limits.Y,'YTick',objP.Axes.Tick.YTick,'YTickLabel',objP.Axes.Tick.YTickLabel,'TickLength',[0,0]); set(get(gca,'YAxis'),'Fontsize',objP.Axes.FontSize{1}(2)); ylabel(objP.Axes.Label(2),'Fontsize',objP.Axes.FontSize{2}(2)); ytickangle(0);
            pbaspect(objP.pbaspect);
            %  > Interpreter.
            fn = fieldnames(get(groot,'factory'));
            fc = find(contains(fn,'Interpreter'));
            for i = 1:numel(fc)
                set(groot,strrep(fn{fc(i)},'factory','Default'),'Latex');
            end
        end
        %  > 1.3.2. -------------------------------------------------------
        %  > Auxiliary function #2.
        function [] = AF_2(objP,Z)
            %  > Axes.
            str.X = ["XTick","XTickLabel"];
            str.Y = ["YTick","YTickLabel"];
            str.Z = ["ZTick","ZTickLabel"];
            for i = ["X","Y","Z"]
                if isempty(objP.Axes.Tick.(str.(i)(2)))
                    if i == "X" || i == "Y"
                        objP.Axes.Tick.(str.(i)(2)) = cell(numel(objP.Axes.Tick.(str.(i)(1))),1);
                        for j = 1:size(objP.Axes.Tick.(str.(i)(2)),1)
                            objP.Axes.Tick.(str.(i)(2)){j} = join(["$",num2str(objP.Axes.Tick.(str.(i)(1))(j)),"$"],'');
                        end
                    else
                        objP.Axes.Limits.Z = [min(Z,[],1),max(Z,[],1)];
                        if  objP.Axes.Limits.Z(1) < objP.ET
                            objP.Axes.Limits.Z(1) = objP.ET;
                        else
                            objP.Axes.Limits.Z(1) = 10.^(ceil(log10(objP.Axes.Limits.Z(1)))-1);
                        end
                        if  objP.Axes.Limits.Z(2) < objP.ET
                            objP.Axes.Limits.Z(2) = objP.ET;
                        else
                            objP.Axes.Limits.Z(2) = 10.^(ceil(log10(objP.Axes.Limits.Z(2)))+0);
                        end
                        if  objP.Axes.Limits.Z(1) == objP.Axes.Limits.Z(2)
                            objP.Axes.Limits.Z(1) = objP.Axes.Limits.Z(2)./10;
                        end
                        if ~isempty(objP.Axes.Tick.(str.(i)(1)))
                            objP.Axes.Tick.(str.(i)(2)) = cell(numel(objP.Axes.Tick.(str.(i)(1))),1);
                            for j = 1:size(objP.Axes.Tick.(str.(i)(2)),1)
                                objP.Axes.Tick.(str.(i)(2)){j} = join(["$",num2str(objP.Axes.Tick.(str.(i)(1))(j)),"$"],'');
                            end
                        else
                            ZLim                        = log10(objP.Axes.Limits.(i));
                            objP.Axes.Tick.(str.(i)(2)) = cell (diff(ZLim)+1,1);
                            for j = 1:size(objP.Axes.Tick.(str.(i)(2)),1)
                                if ZLim(1)-1+j < 0
                                    if numel(num2str(abs(ZLim(1)-1+j))) == 1
                                        objP.Axes.Tick.(str.(i)(2)){j} = join(["$10^{-0",num2str(abs(ZLim(1)-1+j)),"}$"],'');
                                    else
                                        objP.Axes.Tick.(str.(i)(2)){j} = join(["$10^{- ",num2str(abs(ZLim(1)-1+j)),"}$"],'');
                                    end
                                else
                                    if numel(num2str(abs(ZLim(1)-1+j))) == 1
                                        objP.Axes.Tick.(str.(i)(2)){j} = join(["$10^{+0",num2str(abs(ZLim(1)-1+j)),"}$"],'');
                                    else
                                        objP.Axes.Tick.(str.(i)(2)){j} = join(["$10^{+ ",num2str(abs(ZLim(1)-1+j)),"}$"],'');
                                    end
                                end
                                objP.Axes.Tick.(str.(i)(1))(j) = 10.^(ZLim(1)-1+j);
                            end
                            if ~objP.Axes.Tick.(i)
                                nrmv.(i) =  true(size(objP.Axes.Tick.(str.(i)(1))));
                            else
                                nrmv.(i) = ~mod(log10(objP.Axes.Tick.(str.(i)(1))),2);
                            end
                            objP.Axes.Tick.(str.(i)(1)) = objP.Axes.Tick.(str.(i)(1))(nrmv.(i));
                            objP.Axes.Tick.(str.(i)(2)) = objP.Axes.Tick.(str.(i)(2))(nrmv.(i));
                        end
                    end
                end
            end
            clim(gca,objP.Axes.Limits.Z);
            if objP.Colorbar.Scale
                set(gca,'ColorScale','Log');
            end
            set  (gca,'Box','On','Clipping','On');
            set  (gca,'XLim',objP.Axes.Limits.X,'XTick',objP.Axes.Tick.XTick,'XTickLabel',objP.Axes.Tick.XTickLabel,'TickLength',[0,0]); set(get(gca,'XAxis'),'Fontsize',objP.Axes.FontSize{1}(1)); xlabel(objP.Axes.Label(1),'Fontsize',objP.Axes.FontSize{2}(1)); xtickangle(0);
            set  (gca,'YLim',objP.Axes.Limits.Y,'YTick',objP.Axes.Tick.YTick,'YTickLabel',objP.Axes.Tick.YTickLabel,'TickLength',[0,0]); set(get(gca,'YAxis'),'Fontsize',objP.Axes.FontSize{1}(2)); ylabel(objP.Axes.Label(2),'Fontsize',objP.Axes.FontSize{2}(2)); ytickangle(0);
            pbaspect(objP.pbaspect);
            %  > Colorbar.
            if objP.Colorbar.Display
                c  = colorbar(gca,'Box','On','Fontsize',objP.Colorbar.FontSize,'Ticks',objP.Axes.Tick.ZTick,'YTickLabel',objP.Axes.Tick.ZTickLabel);
                bb = get(c,'Position');
                     set(c,'Position',[bb(1),bb(2),0,bb(4)]+objP.Colorbar.Position);
            end
            if objP.Colorbar.Map
                AdvancedColormap('Thermal');
            else
                colormap(objP.C(unique(Z,'rows'),:));
            end
            %  > Interpreter.
            fn = fieldnames(get(groot,'factory'));
            fc = find(contains(fn,'Interpreter'));
            for i = 1:numel(fc)
                set(groot,strrep(fn{fc(i)},'factory','Default'),'Latex');
            end
        end
        % >> 1.4. ---------------------------------------------------------
        %  > 1.4.1. -------------------------------------------------------
        %  > Auxiliary function #1.
        function [Lf] = lf(c,n,dX)
            mu      = cat(1,mean(dX,1),std(dX,0,1));
            dXn     = (dX-mu(1,:))./mu(2,:);
            Df      = c(1).c(1,:).*dXn(:,1).^c(1).e(1,:).*dXn(:,2).^c(1).e(2,:);
            DF      = c(1).c(1,:).*dX (:,1).^c(1).e(1,:).*dX (:,2).^c(1).e(2,:);
            sqrtWf  = sqrt(sqrt(sum(dXn.^2,2)).^(-n));
            sqrtWDf = sqrtWf.*Df;
            mf      = B1_2D.PR2(c);
            Mf      = B1_2D.StM(c,mf,mu,{(sqrtWf.*matlab.internal.math.leastSquaresFit(sqrtWDf,eye(size(sqrtWDf,1)))')'});
            Pf      = Mf{1};
            Lf      = DF*Pf;
        end
    end
end