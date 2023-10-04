classdef B1_2D
    methods (Static)
        %% > 1. -----------------------------------------------------------
        % >> 1.1. ---------------------------------------------------------
        %  > Initialize field 's'.
        function [s] = Initialize_s(inp,msh)
            %  > Auxiliary variables.
            a = msh.c.Nc;            %  > Nc.
            b = msh.f.Nf;            %  > Nf.
            c = size( msh.c.f.if,2); %  > Nf/cell.
            d = sum (~msh.f.l,1);
            if ~inp.PS.EE
                fn = ["c"];
            else
                fn = ["c","r"];
            end

            %  > Loop through...
            for i = fn
                %  > Field: 's'.
                switch i
                    case "c", e = 1;
                    case "r", e = 2;
                end
                s.(i).s.M       = cell (b,e);
                s.(i).s.Q.c     = cell (b,e);
                s.(i).s.Q.m     = cell (b,e);
                s.(i).s.Q.q.x   = cell (b,e);
                s.(i).s.Q.q.w   = cell (b,e);
                s.(i).s.Q.q.v   = cell (b,e);
                if inp.T == 2 && i == "c"
                    for k = 1:msh.f.Nf
                        s.(i).s.sp.c{k,1} = reshape  (msh.f.ic{k,1},[],1);
                        s.(i).s.sp.f{k,1} = RunLength(sort(reshape(msh.c.f.if(msh.f.ic{k,1},:),[],1),1));
                        s.(i).s.ss.c{k,1} = AuxiliaryFunctions.setdiff(RunLength(sort(reshape(cat(2,msh.f.ic{s.(i).s.sp.f{k,1},1}),[],1),1)),s.(i).s.sp.c{k,1});
                        s.(i).s.ss.f{k,1} = AuxiliaryFunctions.setdiff(RunLength(sort(reshape(    msh.c.f.if(s.(i).s.ss.c{k,1},:) ,[],1),1)),s.(i).s.sp.f{k,1});
                    end
                end
                s.(i).s.st      = cell (b,e);
                %  > Field: 'm'.
                s.(i).m.A       = zeros(a);
                s.(i).m.b       = zeros(a,1);
                s.(i).m.nnz.At  = 0;
                %  > Field: 'x.s'.
                s.(i).x.s.nv.a  = zeros(a,1);
                s.(i).x.s.nv.x  = zeros(a,1);
                s.(i).x.s.tfunc = cell (b,e);
                s.(i).x.s.wfunc = cell (d,e);
                if i == "c"
                    %  > Field: 'x.x'.
                    if ~inp.PS.EE
                        for j = ["a","x"]
                            switch j
                                case "a"
                                    s.(i).x.x.dx   .(j)    = zeros(a,c,3); % >> w/ analytical solution.
                                    s.(i).x.x.vf   .(j)    = cell (b,1);   %  > A /LO.
                                    s.(i).x.x.xcunc.(j)    = zeros(a,c,1); %  > A /LO.
                                    s.(i).x.x.xfunc.(j)    = zeros(b,2,1); %  > A /LO.
                                case "x"
                                    s.(i).x.x.vf   .(j)    = cell (b,1);   %  > LO/LO.
                                    s.(i).x.x.xcunc.(j)    = zeros(a,c,1); %  > LO/LO.
                                    s.(i).x.x.xfunc.(j)    = zeros(b,2,1); %  > LO/LO.
                            end
                        end
                    else
                        for j = ["a","x"]
                            switch j
                                case "a"
                                    s.(i).x.x.dx   .(j)    = zeros(a,c,1); % >> w/ analytical solution.
                                    s.(i).x.x.vf   .(j)    = cell (b,1);   %  > A /LO.
                                    s.(i).x.x.xcunc.(j)    = zeros(a,c,1); %  > A /LO.
                                    s.(i).x.x.xfunc.(j)    = zeros(b,2,1); %  > A /LO.
                                case "x"
                                    s.(i).x.x.dx   .(j)    = cell (1,2);   % >> w/ discrete solution.
                                    s.(i).x.x.dx   .(j){1} = zeros(a,c,1); %  > w/ discrete solution (for tau-extrapolation...).
                                    s.(i).x.x.dx   .(j){2} = zeros(a,c,3); %  > w/ discrete solution (for HO estimate...).
                                    s.(i).x.x.vf   .(j)    = cell (b,4);   %  > LO/LO,LO/HO,HO/LO,HO/HO.
                                    s.(i).x.x.xcunc.(j)    = zeros(a,c,4); %  > LO/LO,LO/HO,HO/LO,HO/HO.
                                    s.(i).x.x.xfunc.(j)    = zeros(b,2,4); %  > LO/LO,LO/HO,HO/LO,HO/HO.
                            end
                        end
                    end
                end
            end
            %  > Field: 'u'.
            for i = fn
                switch i
                    case "c"
                        s.u.p.(i){1} = repmat(repmat(inp.PS.P,1,size(inp.C,2))  ,b,1);
                        s.u.s.(i){1} = (1:b)';
                        if inp.T == 2
                            s.u.r = [];
                        end
                    case "r"
                        s.u.p.(i){1} = repmat(repmat(inp.PS.P,1,size(inp.C,2))+2,b,1);
                        s.u.p.(i){2} = repmat(repmat(inp.PS.P,1,size(inp.C,2))+2,b,1);
                        s.u.s.(i){1} = (1:b)';
                        s.u.s.(i){2} = (1:b)';
                end
            end
        end
        % >> 1.2. ---------------------------------------------------------
        %  > Update field 's'.
        function [s] = Update_s(inp,msh,f,s)
            s = B1_2D.Update_ss(inp,msh,f,s);
            s = B1_2D.Update_sm(    msh,f,s);
            s = B1_2D.Update_sx(inp,msh,f,s);
        end
        % >> 1.3. ---------------------------------------------------------
        %  > Discretization stencil setup/flux reconstruction.
        function [s] = Update_ss(inp,msh,f,s)
            %  > Loop through...
            for i = reshape(setdiff(string(fieldnames(s)),"u"),1,[])
                %  > Initialize...
                for j = 1:size(s.u.p.(i),2)
                    pj.(i){j,1} = AuxiliaryFunctions.unique_r(s.u.p.(i){j}(s.u.s.(i){j},:));
                end
                p.(i) = AuxiliaryFunctions.unique_r(cat(1,pj.(i){:,1}));
                for j = 1:size(p.(i),1)
                    Qc.(i){j,1} = B1_2D.PR1 ( p.(i)(j,:)); %  > Polynomial regression coefficients.
                    Qq.(i){j,1} = A3_2D.Q_1D( p.(i)(j,:)); %  > 1D Quadrature.
                end
                for j = 1:size(p.(i),1)
                    Qm.(i){j,1} = B1_2D.PR2 (Qc.(i){j,1}); %  > Nchoosek matrix.
                end
                %  > Update polynomial regression coefficients...
                for j = 1:size(s.u.p.(i),2)
                    for k = s.u.s.(i){j}'
                        s.(i).s.Q.c{k,j} = Qc.(i){all(s.u.p.(i){j}(k,:) == p.(i),2)};
                        s.(i).s.Q.m{k,j} = Qm.(i){all(s.u.p.(i){j}(k,:) == p.(i),2)};
                    end
                end
                %  > Update 1D Quadrature...
                for j = 1:size(s.u.p.(i),2)
                    for k = s.u.s.(i){j}'
                        s.(i).s.Q.q.x{k,j} = f.qd(Qq.(i){all(s.u.p.(i){j}(k,:) == p.(i),2)}.x',msh.f.xy.v{k,1});
                        s.(i).s.Q.q.w{k,j} =      Qq.(i){all(s.u.p.(i){j}(k,:) == p.(i),2)}.w;
                        for l = 1:size(inp.C,1)
                            for m = 1:size(inp.C,2)
                                s.(i).s.Q.q.v{k,j}{l,m} = inp.C{l,m}(s.(i).s.Q.q.x{k,j});
                            end
                        end
                    end
                end
                
                % >> Discretization stencil setup.
                %  > s.
                [s.(i),s.u.p.(i)] = B1_2D.Ds(inp,msh,f,s.(i),s.u.p.(i),s.u.s.(i));
                % >> Flux reconstruction.
                if ~any(any(cellfun(@isempty,s.(i).s.M),1),2)
                    for j = 1:size(s.(i).x.s.tfunc,2)
                        for k = s.u.s.(i){j}'
                            %  > C.
                            for l = ["c","e"]
                                C{1}.(l){1} = s.(i).s.Q.c{k,j}(1).(l);
                                C{1}.(l){2} = s.(i).s.Q.c{k,j}(1).(l);
                            end
                            C{2} = s.(i).s.Q.c{k,j}(2);
                            %  > wVdfunc.
                            wVdfunc = cell(1,size(inp.C,2));
                            for l = 1:size(wVdfunc,2)
                                for m = 1:size(wVdfunc,2)
                                    wVdfunc{1,l}(m,:) = s.(i).s.Q.q.w{k,j}*[s.(i).s.Q.q.v{k,j}{l,m}.*C{l}.c{m}.*(s.(i).s.Q.q.x{k,j}(:,1)-msh.f.xy.c(k,1)).^C{l}.e{m}(1,:).*(s.(i).s.Q.q.x{k,j}(:,2)-msh.f.xy.c(k,2)).^C{l}.e{m}(2,:)];
                                end
                            end
                            wVdfuncT = sum(cat(3,wVdfunc{1,:}),3);
                            %  > wfunc (if required) and tfunc.
                            if ~msh.f.l(k,1)
                                s.(i).x.s.wfunc{f.bd.i == k,j} = wVdfuncT*s.(i).s.M{k,j}.kf;
                            end
                            s.(i).x.s.tfunc{k,j} = wVdfuncT*s.(i).s.M{k,j}.Pf;
                        end
                    end
                end
            end
        end
        %  > 1.3.1. -------------------------------------------------------
        %  > Discretization stencil setup.
        function [ss,su_p] = Ds(inp,msh,f,ss,su_p,su_s)
            for i = 1:size(ss.x.s.tfunc,2)
                if i == 1 && size(ss.x.s.tfunc,2) ~= 1
                    dp = false;
                else
                    dp = true;
                end
                for j = su_s{i}'; p = su_p{i}(j,:);
                                  q = sum(~msh.f.l(j,1),1).*max(ceil(p./2),[],2);
                    %  > Initialize...
                    bf     = 0;
                    Cf     = 0;
                    d      = false(size(p));
                    F.Cond = false;
                    F.Stop = false;
                    xf     = msh.f.xy.c(j,:);

                    %  > Assemble matrix bf (if required)...
                    if ~msh.f.l(j,1)
                        bf = B1_2D.Assemble_bf(f.fh,f.bd.Nf(f.bd.i == j,:),f.bd.t(f.bd.i == j,1),ss.s.Q.q.v{j,i},ss.s.Q.q.x{j,i});
                    end
                    %  > Initialize...
                    msh_d = {msh.d{1}(j,:),msh.d{2}};
                    if i == 1
                        for k = ["i","X"]
                            st.(k) = [];
                        end
                        st.l = false(0);
                    else
                        if ~isempty(ss.s.st{j,1})
                            for k = ["i","l","X"]
                                st.(k) = ss.s.st{j,1}.(k);
                            end
                            msh_d{1} = msh_d{1}(1,~(any(msh_d{2}(1,:) == st.i & msh_d{2}(2,:) == int8(st.l),1)));
                            msh_d{2} = msh_d{2}(:,~(any(msh_d{2}(1,:) == st.i & msh_d{2}(2,:) == int8(st.l),1)));
                        else
                            continue;
                        end
                    end
                    if dp
                        if ~msh.f.l(j,1)
                            st_i = [st.i(~st.l);j];
                        else
                            st_i = [st.i(~st.l)];
                        end
                        if ~isempty(st_i)
                            [d,u] = B1_2D.fd_2(msh,d,st_i);
                            if any(u,2)
                                p       (u)   = p(u)+1;
                                ss.s.Q.c{j,i} = B1_2D.PR1(p);
                            end
                        end
                    end
                    %  > Loop through...
                    while 1
                        while 1
                            if ~F.Cond && size(st.i,1)+q >= size(ss.s.Q.c{j,i}(1).e,2)
                                break;
                            else
                                ns = size(ss.s.Q.c{j,i}(1).e,2)-(size(st.i,1)+q);
                                if ns <= 0
                                    ns = 1;
                                end
                                [msh_d,F.Stop,si,sl,xy] = ...
                                    B1_2D.Add_ds(msh_d,ns,msh.c.c.xy.c,msh.f.xy.c);
                                [st.i] = [st.i;si];
                                [st.l] = [st.l;sl];
                                [st.X] = [st.X;xy];
                                if ~F.Stop
                                    if dp && ~all(d,2) && any(~sl,1)
                                        [d,u] = B1_2D.fd_2(msh,d,si(~sl));
                                        if any(u,2)
                                            p       (u)   = p(u)+1;
                                            ss.s.Q.c{j,i} = B1_2D.PR1(p);
                                        end
                                    end
                                    if size(st.i,1)+q >= size(ss.s.Q.c{j,i}(1).e,2)
                                        break;
                                    end
                                else
                                    break;
                                end
                            end
                        end
                        if ~F.Stop
                            %  > Assemble matrices Cf (if required) and sqrtWDf...
                            dx      = st.X-xf;
                            mu      = cat (1,mean(dx,1),std(dx,0,1));
                            sqrtWf  = sqrt(inp.F.Wf(dx./mu(2,:),2));
                            sqrtWDf = sqrtWf.*B1_2D.Assemble_Df(inp,ss.s.Q.c{j,i},f,st.i,st.l,(dx-mu(1,:))./mu(2,:),st.X);
                            if ~msh.f.l(j,1)
                                Cf = B1_2D.Assemble_Cf(ss.s.Q.c{j,i},f.bd.Nf(f.bd.i == j,:),f.bd.t(f.bd.i == j,1),(ss.s.Q.q.x{j,i}-xf-mu(1,:))./mu(2,:),ss.s.Q.q.v{j,i});
                            end
                            %  > Check matrix (column) rank...
                            if ~msh.f.l(j,1)
                                X = [sqrtWDf;Cf];
                            else
                                X = [sqrtWDf];
                            end
                            if rank(X) < size(ss.s.Q.c{j,i}(1).e,2)
                                if (size(sqrtWDf,1)-q)/size (ss.s.Q.c{j,i}(1).e,2) > 3.0
                                    F.Cond = false;
                                    F.Stop = true;
                                else
                                    F.Cond = true;
                                end
                            else
                                %  > Update (if required)...
                                if any(su_p{i}(j,:) ~= p,2)
                                    m = all(su_p{i} == p,2);
                                    if all(~m,1)
                                        ss.s.Q.m{j,i} = B1_2D.PR2(ss.s.Q.c{j,i});
                                    else
                                        ss.s.Q.m{j,i} = ss.s.Q.m{find(m,1),i};
                                    end
                                end
                                %  > Assemble matrices kf (if required) and Pf...
                                if ~msh.f.l(j,1)
                                    [A,B]               = qr([sqrtWDf;Cf],0);
                                    [L]                 =    [true(size(sqrtWDf,1),1);false(size(Cf,1),1)];
                                    [M]                 = B1_2D.StM(ss.s.Q.c{j,i},ss.s.Q.m{j,i},mu,{matlab.internal.math.leastSquaresFit(B,(matlab.internal.math.leastSquaresFit(A(~L,:)',eye(size(A(~L,:)',1))))'*bf),...
                                                                                                    matlab.internal.math.leastSquaresFit(B,A(L,:)'-A(~L,:)'*matlab.internal.math.leastSquaresFit(A(~L,:)',eye(size(A(~L,:)',1)))*A(L,:)').*sqrtWf'});
                                    ss.s.M{j,i}.bf      = bf;
                                    ss.s.M{j,i}.Cf      = Cf;
                                    ss.s.M{j,i}.sqrtWDf = sqrtWDf;
                                    ss.s.M{j,i}.kf      = M{1};
                                    ss.s.M{j,i}.Pf      = M{2};
                                else
                                    [M]                 = B1_2D.StM(ss.s.Q.c{j,i},ss.s.Q.m{j,i},mu,{(sqrtWf.*matlab.internal.math.leastSquaresFit(sqrtWDf,eye(size(sqrtWDf,1)))')'});
                                    ss.s.M{j,i}.sqrtWDf = sqrtWDf;
                                    ss.s.M{j,i}.Pf      = M{1};
                                end
                                ss.s.st{j,i}    = st;
                                su_p   {i}(j,:) = p;
                                break;
                            end
                        else
                            ss.s.M {j,i} = [];
                            ss.s.st{j,i} = [];
                            break;
                        end
                    end
                end
            end
        end
        %  > 1.3.1.1. -----------------------------------------------------
        %  > Auxiliary function #1 (local  direction).
        function [d] = fd_1(N)
            d = [abs(N(:,1)) >= abs(N(:,2)),...
                 abs(N(:,2)) >= abs(N(:,1))];
        end
        %  > 1.3.1.2. -----------------------------------------------------
        %  > Auxiliary function #2 (global direction).
        function [d,u] = fd_2(msh,d,sf)
            Nf = zeros(size(sf,1),size(d,2));
            for i = 1:size(Nf,1)
                Nf(i,:) = msh.c.f.Nf{msh.f.ic{sf(i),1},1}(msh.c.f.if(msh.f.ic{sf(i,1),1},:) == sf(i,1),:);
            end
            u = ~d & any(B1_2D.fd_1(Nf),1);
            d =  d | u;
        end
        %  > 1.3.1.3. -----------------------------------------------------
        %  > Assemble matrix Df.
        function [Df] = Assemble_Df(inp,c,f,si,sl,dx,xt)
            %  > Check boundary type...
            t        = f.bd.t( AuxiliaryFunctions.find_c(f.bd.i,si(~sl,1)),1);
            bf       = find  (~sl);
            %  > ...and assemble matrix Df accordingly.
            Df(sl,:) = c(1).c(1,:).*dx(sl,1).^c(1).e(1,:).*dx(sl,2).^c(1).e(2,:);
            if ~isempty(t)
                for i = bf'
                    switch t(bf == i)
                        case "Dirichlet"
                            Df(i,:) = c(1).c(1,:).*dx(i,1).^c(1).e(1,:).*dx(i,2).^c(1).e(2,:);
                        case {"Neumann","Robin"}
                            n  = f.bd.Nf(f.bd.i == si(i,1),:);
                            D2 = zeros  (size(n,2),size(c(1).c,2));
                            switch t(bf == i)
                                case "Neumann"
                                    D2 = zeros(size(n,2),size(c(2).c{j},2));
                                    for j = 1:size(n,2)
                                        D2(j,:) = c(2).c{j}(1,:).*dx(i,1).^c(2).e{j}(1,:).*dx(i,2).^c(2).e{j}(2,:);
                                    end
                                case "Robin"
                                    D1 = repmat(c(1).c(1,:).*dx(i,1).^c(1).e(1,:).*dx(i,2).^c(1).e(2,:),size(n,2),1);
                                    for j = 1:size(n,2)
                                        D2(j,:) = inp.C{2,j}(xt(i,:))./inp.C{1,j}(xt(i,:)).*c(2).c{j}(1,:).*dx(i,1).^c(2).e{j}(1,:).*dx(i,2).^c(2).e{j}(2,:);
                                    end
                                    D2(isinf(D2) | isnan(D2)) = 0; D2 = D1+D2;
                            end
                            Df(i,:) = n*D2;
                    end
                end
            end
        end
        %  > 1.3.1.4. -----------------------------------------------------
        %  > Assemble (or add to...) discretization stencil.
        function [msh_d,F,si,sl,xy] = Add_ds(msh_d,m,xy_c,xy_f)
            n = m;
            if(any(msh_d{1} == 0,2))
                n = m+1;
            end
            dm = mink(msh_d{1},n);
            e  = msh_d{1} <= dm(end) & msh_d{1} ~= 0;
            F  = false;
            si = msh_d{2}(1,e)';
            sl = msh_d{2}(2,e)'; sl = logical(sl);
            xy = [xy_c(si( sl),:);
                  xy_f(si(~sl),:)];
            if isempty(si)
                F = true;
            end
            msh_d{1} = msh_d{1}(1,~e);
            msh_d{2} = msh_d{2}(:,~e);
        end
        %  > 1.3.1.5. -----------------------------------------------------
        %  > Assemble matrix bf.
        function [bf] = Assemble_bf(fh,n,t,vg,xg)
            switch t
                case "Dirichlet"
                    bf = fh.f(xg);
                case {"Neumann","Robin"}
                    V2 = zeros(size(xg));
                    switch t
                        case "Neumann"
                            for i = 1:size(n,2)
                                V2(:,i) = fh.d{i}(xg);
                            end
                        case "Robin"
                            V1 = repmat(fh.f(xg),1,size(n,2));
                            for i = 1:size(n,2)
                                V2(:,i) = vg{2,i}./vg{1,i}.*fh.d{i}(xg);
                            end
                            V2(isinf(V2) | isnan(V2)) = 0; V2 = V1+V2;
                    end
                    bf = V2*n';
            end
        end
        %  > 1.3.1.6. -----------------------------------------------------
        %  > Assemble matrix Cf.
        function [Cf] = Assemble_Cf(c,n,t,dxg,vg)
            switch t
                case "Dirichlet"
                    Cf = c(1).c(1,:).*dxg(:,1).^c(1).e(1,:).*dxg(:,2).^c(1).e(2,:);
                case {"Neumann","Robin"}
                    D2 = cell(size(n));
                    switch t
                        case "Neumann"
                            for i = 1:size(n,2)
                                D2{i} = n(1,i).*c(2).c{i}(1,:).*dxg(:,1).^c(2).e{i}(1,:).*dxg(:,2).^c(2).e{i}(2,:);
                            end
                        case "Robin"
                            D1 = c(1).c(1,:).*dxg(:,1).^c(1).e(1,:).*dxg(:,2).^c(1).e(2,:);
                            for i = 1:size(n,2)
                                D2{i} = vg{2,i}./vg{1,i}.*c(2).c{i}(1,:).*dxg(:,1).^c(2).e{i}(1,:).*dxg(:,2).^c(2).e{i}(2,:); D2{i}(isinf(D2{i}) | isnan(D2{i})) = 0;
                                D2{i} = n(1,i).*(D1+D2{i});
                            end
                    end
                    Cf = sum(cat(3,D2{:}),3);
            end
        end
        %  > 1.3.2. -------------------------------------------------------
        %  > Polyfit standardization.
        function [M3] = StM(c,t,mu,M1)
            %  > Initialize...
            M2 = zeros(size(t.l));
            M3 = cell (size(M1));
            for i = 1:size(M3,2)
                M3{i} = zeros(size(M1{i}));
            end
            
            for i = 1:size(c(1).e,2)
                j              = c(1).e(:,t.l(i,:));
                k              = c(1).e(:,t.l(i,:))-c(1).e(:,i);
                M2(i,t.l(i,:)) = ...
                    t.m(i,t.l(i,:)).*(-mu(1,1)).^k(1,:)./(mu(2,1).^j(1,:)).*(-mu(1,2)).^k(2,:)./(mu(2,2).^j(2,:));
            end
            for l = 1:size(M3,2)
                M3{l} = M2*M1{l};
            end
        end
        % >> 1.4. ---------------------------------------------------------
        %  > Compute polynomial regression coefficients/exponents.
        %  > 1.4.1. -------------------------------------------------------
        function [t] = PR1(p)
            %  > Auxiliary variables.
            m      = max(p,[],2);
            n{1,2} = ones(m+1,1)*(0:m);
            n{1,1} = n{2}';
            o{2,2} = ones(m+1,1)*[0,0:m-1];
            o{1,1} = o{2,2}';
            o{1,2} = ones(m+1,1)*(0:m);
            o{2,1} = o{1,2}';
            c{1}   = o{2,1};
            c{2}   = o{1,2};
            s      = sum(cat(3,n{:}),3);
            l      = n{1} <= p(1,1) & n{2} <= p(1,2) & s <= m;
            if p(1) ~= p(2)
                l(s == m) = n{p == m}(s == m) >= n{p ~= m}(s == m); 
            end
            
            %  > t(1).
            t(1).c      = ones(1,sum(l,"all"));  %  > c.
            t(1).e(1,:) = n{1}(l);               %  > e(x).
            t(1).e(2,:) = n{2}(l);               %  > e(y).
            %  > t(2).
            for k = 1:numel(c)
                t(2).c{k}(1,:) = c{k}  (l);      %  > c.
                t(2).e{k}(1,:) = o{k,1}(l);      %  > e(x).
                t(2).e{k}(2,:) = o{k,2}(l);      %  > e(y).
                u              = t(2).c{k} == 0;
                t(2).e{k}(:,u) = 0;
            end
        end
        %  > 1.4.2. -------------------------------------------------------
        function [t] = PR2(c)
            %  > Initialize...
            x   = 1:size(c(1).e,2);
            t.l = repmat(c(1).e(1,:),size(c(1).e,2),1) >= c(1).e(1,:)' & repmat(c(1).e(2,:),size(c(1).e,2),1) >= c(1).e(2,:)';
            t.m = zeros (size(c(1).e,2));
            
            for i = 1:size(c(1).e,2)
                j = c(1).e(:,t.l(i,:));
                k = c(1).e(:,t.l(i,:))-c(1).e(:,i);
                l = x(t.l(i,:));
                for m = 1:sum(t.l(i,:),2)
                    t.m(i,l(m)) = prod(((j(1,m)-k(1,m)+1):j(1,m))./(1:k(1,m))).*prod(((j(2,m)-k(2,m)+1):j(2,m))./(1:k(2,m)));
                end
            end
        end
        % >> 1.5. ---------------------------------------------------------
        %  > Update matrix A, vector b and nodal solution.
        function [s] = Update_sm(msh,f,s)
            for i = reshape(setdiff(string(fieldnames(s)),"u"),1,[])
                if ~any(any(cellfun(@isempty,s.(i).s.M),1),2)
                    %  > Select...
                    switch i
                        case "c", x = 1;
                        case "r", x = 2;
                    end
                    %  > Reset if...
                    r = RunLength(sort(cat(2,msh.f.ic{s.u.s.(i){x},1}),2));
                    if any(~s.(i).m.A,"all")
                        s.(i).m.A(r,:) = 0;
                        s.(i).m.b(r,1) = sum(f.st(r,:),2);
                    end
                    %  > For each cell...
                    for j = r
                        %  > For each face...
                        for k = 1:size(msh.c.f.if,2), l = msh.c.f.if(j,k);
                            %  > Auxiliary variables.
                            m  = s.(i).s.st{l,x}.l;
                            n  = s.(i).s.st{l,x}.i( m);
                            o  = s.(i).s.st{l,x}.i(~m);
                            Sf = msh.c.f.Sf{j}    ( k,:);
                            %  > A.
                            s.(i).m.A(j,n) = s.(i).m.A(j,n)+Sf*s.(i).x.s.tfunc{l,x}(:,m);
                            %  > b.
                            if any(~m)
                                s.(i).m.b(j,1) = s.(i).m.b(j,1)-Sf*s.(i).x.s.tfunc{l,x}(:,~m)*f.bd.v(AuxiliaryFunctions.find_c(f.bd.i,o),1);
                            end
                            %  > Add to RHS...
                            if ~msh.f.l(l,1)
                                s.(i).m.b(j,1) = s.(i).m.b(j,1)-Sf*s.(i).x.s.wfunc{f.bd.i == l,x};
                            end
                        end
                    end
                    %  > nnz.
                    s.(i).m.nnz.At      = nnz(s.(i).m.A);
                    %  > nv.
                    s.(i).x.s.nv.a(:,1) = f.fh.f(msh.c.c.xy.c);
                    s.(i).x.s.nv.x(:,1) = s.(i).m.A\s.(i).m.b;
                end
            end
        end
        % >> 1.6. ---------------------------------------------------------
        %  > Update facial values, etc.
        function [s] = Update_sx(inp,msh,f,s)
            if ~any(any(cellfun(@isempty,s.c.s.M),1),2)
                %  > Auxiliary variables.
                iu.a = s.u.s.c{1}';
                iu.x = 1:msh.f.Nf;
                sl   = s.c;
                nv.c = s.c.x.s.nv;
                
                %  > vf.
                s.c.x.x.vf.a   (iu.a,1)   = B1_2D.vf_i(    f,iu.a,sl,nv.c.a);                       %  > A /LO.
                s.c.x.x.vf.x   (iu.x,1)   = B1_2D.vf_i(    f,iu.x,sl,nv.c.x);                       %  > LO/LO.
                %  > xfunc.
                s.c.x.x.xfunc.a(iu.a,:,1) = B1_2D.xf_i(msh,f,iu.a,sl,s.c.x.x.vf.a(:,1));            %  > A /LO.
                s.c.x.x.xfunc.x(iu.x,:,1) = B1_2D.xf_i(msh,f,iu.x,sl,s.c.x.x.vf.x(:,1));            %  > LO/LO.
                if inp.PS.EE
                    if ~any(any(cellfun(@isempty,s.r.s.M),1),2)
                        sh(size(s.r.s.st,2)) = struct();
                        for i = 1:size(s.r.s.st,2)
                            sh(i).s.M       = s.r.s.M      (:,i);
                            sh(i).s.st      = s.r.s.st     (:,i);
                            sh(i).x.s.tfunc = s.r.x.s.tfunc(:,i);
                            sh(i).x.s.wfunc = s.r.x.s.wfunc(:,i);
                        end
                        nv.r                      = s.r.x.s.nv;
                        s.c.x.x.vf.x   (iu.x,2)   = B1_2D.vf_i(    f,iu.x,sh(1),nv.c.x);            %  > LO/HO (for tau-extrapolation...).
                        s.c.x.x.vf.x   (iu.x,3)   = B1_2D.vf_i(    f,iu.x,sl   ,nv.r.x);            %  > HO/LO (for HO estimate...).
                        s.c.x.x.vf.x   (iu.x,4)   = B1_2D.vf_i(    f,iu.x,sh(2),nv.r.x);            %  > HO/HO.
                        s.c.x.x.xfunc.x(iu.x,:,2) = B1_2D.xf_i(msh,f,iu.x,sh(1),s.c.x.x.vf.x(:,2)); %  > LO/HO (for tau-extrapolation...).
                        s.c.x.x.xfunc.x(iu.x,:,3) = B1_2D.xf_i(msh,f,iu.x,sl   ,s.c.x.x.vf.x(:,3)); %  > HO/LO (for HO estimate...).
                        s.c.x.x.xfunc.x(iu.x,:,4) = B1_2D.xf_i(msh,f,iu.x,sh(2),s.c.x.x.vf.x(:,4)); %  > HO/HO.
                    end
                end
                %  > xcunc.
                for i = 1:msh.f.Nf
                    for j = 1:numel(msh.f.ic{i,1})
                        k = msh.c.f.if(msh.f.ic{i,1}(j),:) == i;
                        for l = ["a","x"]
                            if ismembc(i,iu.(l))
                                for m = 1:size(s.c.x.x.xcunc.(l),3)
                                    s.c.x.x.xcunc.(l)(msh.f.ic{i,1}(j),k,m) = msh.c.f.Sf{msh.f.ic{i,1}(j)}(k,:)*s.c.x.x.xfunc.(l)(i,:,m)';
                                end
                            end
                        end
                    end
                end
                %  > dx (w/ analytical solution)...
                s.c.x.x.dx.a(:,:,1) = s.c.x.x.xcunc.a(:,:,1)-s.c.x.x.xcunc.x(:,:,1);                %  > a0+a1*h^n+...
                s.c.x.x.dx.a(:,:,2) = f.st                  -s.c.x.x.xcunc.x(:,:,1);                %  > a0.
                s.c.x.x.dx.a(:,:,3) = s.c.x.x.xcunc.a(:,:,1)-f.st;                                  %  > a1*h^n+...
                %  > dx (w/ discrete solution)...
                if inp.PS.EE
                    if ~any(any(cellfun(@isempty,s.r.s.M),1),2)
                        %  > ...for tau-extrapolation.
                        s.c.x.x.dx.x{1}(:,:,1) = s.c.x.x.xcunc.x(:,:,1)-s.c.x.x.xcunc.x(:,:,2);     %  > a1*h^n.
                        %  > ...for HO estimate.
                        s.c.x.x.dx.x{2}(:,:,1) = s.c.x.x.xcunc.x(:,:,3)-s.c.x.x.xcunc.x(:,:,1);     %  > a0+a1*h^n.
                        s.c.x.x.dx.x{2}(:,:,2) = s.c.x.x.xcunc.x(:,:,4)-s.c.x.x.xcunc.x(:,:,1);     %  > a0.
                        s.c.x.x.dx.x{2}(:,:,3) = s.c.x.x.xcunc.x(:,:,3)-s.c.x.x.xcunc.x(:,:,4);     %  > a1*h^n.
                    end
                end
            end
       end
        %  > 1.6.1. -------------------------------------------------------
        %  > Auxiliary function #1.
        function [vf] = vf_i(f,iu,s,nv)
            %  > Initialize...
            vf = cell(size(iu,2),1);
            
            %  > For each face...
            for i = 1:size(vf,1), j = iu(1,i); k = s.s.st{j,1}.l;
                %  > Cell value(s).
                vf{i,1}      = zeros(size(k,1),1);
                vf{i,1}(k,1) = nv   (s.s.st{j,1}.i(k,1),1);
                %  > Face value(s).
                if any(~k,1)
                    vf{i,1}(~k,1) = f.bd.v(AuxiliaryFunctions.find_c(f.bd.i,s.s.st{j,1}.i(~k,1)),1);
                end
            end
        end
        %  > 1.6.2. -------------------------------------------------------
        %  > Auxiliary function #2.
        function [xf] = xf_i(msh,f,iu,s,vf)
            %  > Initialize...
            xf = cell(size(iu,2),1);
            
            %  > For each face...
            for i = 1:size(xf,1), j = iu(1,i);
                xf{i,1} = s.x.s.tfunc{j,1}*vf{j,1};
                if ~msh.f.l(j,1)
                    xf{i,1} = xf{i,1}+s.x.s.wfunc{f.bd.i == j,1};
                end
            end
            xf = cat(2,xf{:,1})';
        end
    end
end