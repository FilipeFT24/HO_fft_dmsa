classdef B1_1D
    methods (Static)
        %% > 1. -----------------------------------------------------------
        % >> 1.1. ---------------------------------------------------------
        %  > Initialize field 's'.
        function [s] = Initialize_s(inp,msh)
            %  > Auxiliary variables.
            a = msh.c.Nc;           %  > Nc.
            b = msh.f.Nf;           %  > Nf.
            c = size(msh.c.f.if,2); %  > Nf/cell.
            if ~inp.PS.EE
                fn = ["c"];
            else
                fn = ["c","r"];
            end

            %  > Loop through...
            for i = fn
                %  > Field: 's'.
                switch i
                    case "c", d = 1;
                    case "r", d = 2;
                end
                s.(i).s.M       = cell (b,d);
                if inp.T == 2 && i == "c"
                    for k = 1:msh.f.Nf
                        s.(i).s.sp.c{k,1} = reshape  (msh.f.ic{k,1},[],1);
                        s.(i).s.sp.f{k,1} = RunLength(sort(reshape(msh.c.f.if(msh.f.ic{k,1},:),[],1),1));
                        s.(i).s.ss.c{k,1} = AuxiliaryFunctions.setdiff(RunLength(sort(reshape(cat(2,msh.f.ic{s.(i).s.sp.f{k,1},1}),[],1),1)),s.(i).s.sp.c{k,1});
                        s.(i).s.ss.f{k,1} = AuxiliaryFunctions.setdiff(RunLength(sort(reshape(    msh.c.f.if(s.(i).s.ss.c{k,1},:) ,[],1),1)),s.(i).s.sp.f{k,1});
                    end
                end
                s.(i).s.st      = cell (b,d);
                %  > Field: 'm'.
                s.(i).m.A       = zeros(a);
                s.(i).m.b       = zeros(a,1);
                s.(i).m.nnz.At  = 0;
                %  > Field: 'x.s'.
                s.(i).x.s.nv.a  = zeros(a,1);
                s.(i).x.s.nv.x  = zeros(a,1);
                s.(i).x.s.tfunc = cell (b,d);
                if i == "c"
                    %  > Field: 'x.x'.
                    if ~inp.PS.EE
                        for j = ["a","x"]
                            switch j
                                case "a"
                                    s.(i).x.x.dx   .(j)    = zeros(a,c,3); % >> w/ analytical solution.
                                    s.(i).x.x.vf   .(j)    = cell (b,1);   %  > A /LO.
                                    s.(i).x.x.xcunc.(j)    = zeros(a,c,1); %  > A /LO.
                                    s.(i).x.x.xfunc.(j)    = zeros(b,1);   %  > A /LO.
                                case "x"
                                    s.(i).x.x.vf   .(j)    = cell (b,1);   %  > LO/LO.
                                    s.(i).x.x.xcunc.(j)    = zeros(a,c,1); %  > LO/LO.
                                    s.(i).x.x.xfunc.(j)    = zeros(b,1);   %  > LO/LO.
                            end
                        end
                    else
                        for j = ["a","x"]
                            switch j
                                case "a"
                                    s.(i).x.x.dx   .(j)    = zeros(a,c,3); % >> w/ analytical solution.
                                    s.(i).x.x.vf   .(j)    = cell (b,1);   %  > A /LO,A /HO.
                                    s.(i).x.x.xcunc.(j)    = zeros(a,c,1); %  > A /LO,A /HO.
                                    s.(i).x.x.xfunc.(j)    = zeros(b,1,1); %  > A /LO,A /HO.
                                case "x"
                                    s.(i).x.x.dx   .(j)    = cell (1,2);   % >> w/ discrete solution.
                                    s.(i).x.x.dx   .(j){1} = zeros(a,c,1); %  > w/ discrete solution (for tau-extrapolation...).
                                    s.(i).x.x.dx   .(j){2} = zeros(a,c,3); %  > w/ discrete solution (for HO estimate...).
                                    s.(i).x.x.vf   .(j)    = cell (b,4);   %  > LO/LO,LO/HO,HO/LO,HO/HO.
                                    s.(i).x.x.xcunc.(j)    = zeros(a,c,4); %  > LO/LO,LO/HO,HO/LO,HO/HO.
                                    s.(i).x.x.xfunc.(j)    = zeros(b,4);   %  > LO/LO,LO/HO,HO/LO,HO/HO.
                            end
                        end
                    end
                end
            end
            %  > Field: 'u'.
            for i = fn
                switch i
                    case "c"
                        s.u.p.(i){1} = repmat(inp.PS.P  ,b,1);
                        s.u.s.(i){1} = (1:b)';
                        if inp.T == 2
                            s.u.r = [];
                        end
                    case "r"
                        s.u.p.(i){1} = repmat(inp.PS.P+2,b,1);
                        s.u.p.(i){2} = repmat(inp.PS.P+2,b,1);
                        s.u.s.(i){1} = (1:b)';
                        s.u.s.(i){2} = (1:b)';
                end
            end
        end
        % >> 1.2. ---------------------------------------------------------
        %  > Update field 's'.
        function [s] = Update_s(inp,msh,f,s)
            s = B1_1D.Update_ss(inp,msh,f,s);
            s = B1_1D.Update_sm(    msh,f,s);
            s = B1_1D.Update_sx(inp,msh,f,s);
        end
        % >> 1.3. ---------------------------------------------------------
        %  > Discretization stencil setup/flux reconstruction.
        function [s] = Update_ss(inp,msh,f,s)
            for i = reshape(setdiff(string(fieldnames(s)),"u"),1,[])
                % >> Discretization stencil setup.
                %  > s.
                [s.(i),s.u.p.(i)] = B1_1D.Assemble_Ds(inp,msh,f,s.(i),s.u.p.(i),s.u.s.(i));
                % >> Flux reconstruction.
                if ~any(any(cellfun(@isempty,s.(i).s.M),1),2)
                    for j = 1:size(s.(i).x.s.tfunc,2)
                        for k = s.u.s.(i){j}'
                            %  > Vdfunc.
                            Vdfunc = zeros(size(inp.C,1),size(s.(i).s.M{k,j}.Pf,2));
                            for l = 1:size(Vdfunc,1)
                                Vdfunc(l,l) = inp.C{l}(msh.f.Xv(k,1)).*factorial(l-1);
                            end
                            %  > tfunc.
                            s.(i).x.s.tfunc{k,j} = sum(Vdfunc,1)*s.(i).s.M{k,j}.Pf;
                        end
                    end
                end
            end
        end
        %  > 1.3.1. -------------------------------------------------------
        %  > Discretization stencil setup.
        function [ss,su_p] = Assemble_Ds(inp,msh,f,ss,su_p,su_s)
            x = cell(size(ss.x.s.tfunc,2));
            for i = 1:size(ss.x.s.tfunc,2)
                x{i} = any([(1:msh.f.Nf)' <= ceil(su_p{i}./2),(1:msh.f.Nf)' >= size((1:msh.f.Nf)',1)-(ceil(su_p{i}./2)-1)],2);
                if ~(i == 1 && size(ss.x.s.tfunc,2) ~= 1)
                    for j = su_s{i}', p = su_p{i}(j,1); q = ceil(p./2); sc = cell(1,q);
                        %  > Initialize...
                        if ~x{i}(j)
                            for k = 1:q
                                sc{k} = [j-k,j+k-1];
                            end
                            ss.s.st{j,i}.i = sort(cat(2,sc{:}),2)';
                            ss.s.st{j,i}.l = true(size(ss.s.st{j,i}.i));
                        else
                            l = [j > 1:q;msh.f.Nf-(j-1) > 1:q]';
                            m = any(~l,1);
                            for k = 1:q
                                if l(k,1) && l(k,2)
                                    a = j-k;
                                    b = j+k-1;
                                else
                                    if l(k,2), a = 2.*k-1;
                                    end
                                    if l(k,1), a = msh.f.Nf-2.*k;
                                    end
                                    b = a+1;
                                end
                                sc{k} = [a,b];
                            end
                            if m(1)
                                c              = 1;
                                d              = false;
                                ss.s.st{j,i}.i = cat(2,c,sort(cat(2,sc{:}),2))';
                                ss.s.st{j,i}.l = cat(1,d,true(size(ss.s.st{j,i}.i)-[1,0]));
                            end
                            if m(2)
                                c              = msh.f.Nf;
                                d              = false;
                                ss.s.st{j,i}.i = cat(2,sort(cat(2,sc{:}),2),c)';
                                ss.s.st{j,i}.l = cat(1,true(size(ss.s.st{j,i}.i)-[1,0]),d);
                            end
                        end
                        ss.s.st{j,i}.xt = sort(cat(1,msh.f.Xv(ss.s.st{j,i}.i(~ss.s.st{j,i}.l),1),msh.c.Xc(ss.s.st{j,i}.i(ss.s.st{j,i}.l),1)),1);
                        %  > Assemble matrix Pf...
                        Pf = B1_1D.Assemble_Pf(inp,msh,f,msh.f.Xv(j,1),ss.s.st{j,i}.xt);
                        if ~isempty(Pf)
                            ss.s.M{j,i}.Pf   = Pf;
                            su_p    {i}(j,1) = size(Pf,1)-1;
                        else
                            ss.s.M{j,i}      = []; 
                        end
                    end
                end
            end
            if size(ss.x.s.tfunc,2) ~= 1
                i = 1;
                for j = su_s{i}', p = su_p{i}(j,1); q = ceil(p./2);
                    ss.s.st(j,i) = ss.s.st(j,i+1);
                    if x{i}(j)
                        m = any(~[j > 1:q;msh.f.Nf-(j-1) > 1:q]',1);
                        for k = ["i","l","xt"]
                            ss.s.st{j,i}.(k) = ss.s.st{j,i}.(k);
                            if m(1), ss.s.st{j,i}.(k)(end) = [];
                            end
                            if m(2), ss.s.st{j,i}.(k)(1)   = [];
                            end
                        end
                    end
                    %  > Assemble matrix Pf...
                    Pf = B1_1D.Assemble_Pf(inp,msh,f,msh.f.Xv(j,1),ss.s.st{j,i}.xt);
                    if ~isempty(Pf)
                        ss.s.M{j,i}.Pf   = Pf;
                        su_p    {i}(j,1) = size(Pf,1)-1;
                    else
                        ss.s.M{j,i}      = [];
                    end
                end
            end
        end
        %  > 1.3.2. -------------------------------------------------------
        %  > Assemble matrix Pf.
        function [Pf] = Assemble_Pf(inp,msh,f,xf,xt)
            %  > Auxiliary variables.
            dx = xt-xf;
            p  = 0:size(xt,1)-1;
            l  = ismembc(xt,msh.f.Xv(~msh.f.l,1));
            if any(l,1)
                b = f.bd.i == (find(msh.f.Xv == xt(l,1))); n = f.bd.Nf(b,1);
                                                           t = f.bd.t (b,1);
            end
            
            %  > Polyfit standardization.
            mu = [mean(dx,1),std(dx,0,1)];
            dX = (dx-mu(1))./mu(2);
            %  > Df and Pf.
            if any(l,1) && t ~= "Dirichlet"
                q = p-1; q(q < 0) = 0;
                switch t
                    case "Neumann"
                        Df(~l,:) = dX(~l,1).^p;
                        Df( l,:) = dX( l,1).^q.*p;
                    case "Robin"
                        vg       = inp.C{2}(xf)./inp.C{1}(xf); vg(isinf(vg) | isnan(vg)) = 0;
                        Df(~l,:) = dX(~l,1).^p;
                        Df( l,:) = dX( l,1).^p+vg.*dX(l,1).^q.*p;  
                end
                Df(isinf(Df) | isnan(Df)) = 0; Df(l,:) = n*Df(l,:);
            else
                Df = dX.^p;
            end
            %  > Check matrix (column) rank...
            if rank(Df) < size(Df,2)
                Pf = [];
            else
                Pf = zeros(size(Df));
                PF = eye  (size(Df))/Df;
                for i = 0:p(end)
                    for j = i:p(end)
                        Pf(i+1,:) = Pf(i+1,:)+PF(j+1,:).*nchoosek(j,j-i)./mu(2).^j.*(-mu(1))^(j-i);
                    end
                end
            end
        end
        % >> 1.4. ---------------------------------------------------------
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
                        for k = 1:numel(msh.c.f.if(j,:)), l = msh.c.f.if(j,k);
                            %  > Auxiliary variables.
                            m  = s.(i).s.st{l,x}.l;
                            n  = s.(i).s.st{l,x}.i( m);
                            o  = s.(i).s.st{l,x}.i(~m);
                            Nf = msh.c.f.Nf       ( j,k);
                            %  > A.
                            s.(i).m.A(j,n) = s.(i).m.A(j,n)+Nf.*s.(i).x.s.tfunc{l,x}(m);
                            %  > b.
                            if any(~m)
                                s.(i).m.b(j,1) = s.(i).m.b(j,1)-Nf.*s.(i).x.s.tfunc{l,x}(~m)*f.bd.v(AuxiliaryFunctions.find_c(f.bd.i,o),1);
                            end
                        end
                    end
                    %  > nnz.
                    s.(i).m.nnz.At      = nnz(s.(i).m.A);
                    %  > nv.
                    s.(i).x.s.nv.a(:,1) = f.fh.f(msh.c.Xc);
                    s.(i).x.s.nv.x(:,1) = s.(i).m.A\s.(i).m.b;
                end
            end
        end
        % >> 1.5. ---------------------------------------------------------
        %  > Update facial values, etc.
        function [s] = Update_sx(inp,msh,f,s)
            if ~any(any(cellfun(@isempty,s.c.s.M),1),2)
                %  > Auxiliary variables.
                iu.a = s.u.s.c{1}';
                iu.x = 1:msh.f.Nf;
                sl   = s.c;
                nv.c = s.c.x.s.nv;
                
                %  > vf.
                s.c.x.x.vf.a   (iu.a,1) = B1_1D.vf_i(f,iu.a,sl,nv.c.a);                         %  > A /LO.
                s.c.x.x.vf.x   (iu.x,1) = B1_1D.vf_i(f,iu.x,sl,nv.c.x);                         %  > LO/LO.
                %  > xfunc.
                s.c.x.x.xfunc.a(iu.a,1) = B1_1D.xf_i(  iu.a,sl,s.c.x.x.vf.a(:,1));              %  > A /LO.
                s.c.x.x.xfunc.x(iu.x,1) = B1_1D.xf_i(  iu.x,sl,s.c.x.x.vf.x(:,1));              %  > LO/LO.
                if inp.PS.EE
                    if ~any(any(cellfun(@isempty,s.r.s.M),1),2)
                        sh(size(s.r.s.st,2)) = struct();
                        for i = 1:size(s.r.s.st,2)
                            sh(i).s.st      = s.r.s.st     (:,i);
                            sh(i).x.s.tfunc = s.r.x.s.tfunc(:,i);
                        end
                        nv.r                    = s.r.x.s.nv;
                        s.c.x.x.vf.x   (iu.x,2) = B1_1D.vf_i(f,iu.x,sh(1),nv.c.x);              %  > LO/HO (for tau-extrapolation...).
                        s.c.x.x.vf.x   (iu.x,3) = B1_1D.vf_i(f,iu.x,sl   ,nv.r.x);              %  > HO/LO (for HO estimate...).
                        s.c.x.x.vf.x   (iu.x,4) = B1_1D.vf_i(f,iu.x,sh(2),nv.r.x);              %  > HO/HO.
                        s.c.x.x.xfunc.x(iu.x,2) = B1_1D.xf_i(  iu.x,sh(1),s.c.x.x.vf.x(:,2));   %  > LO/HO (for tau-extrapolation...).
                        s.c.x.x.xfunc.x(iu.x,3) = B1_1D.xf_i(  iu.x,sl   ,s.c.x.x.vf.x(:,3));   %  > HO/LO (for tau-interpolation...).
                        s.c.x.x.xfunc.x(iu.x,4) = B1_1D.xf_i(  iu.x,sh(2),s.c.x.x.vf.x(:,4));   %  > HO/HO.
                    end
                end
                %  > xcunc.
                for i = 1:msh.f.Nf
                    for j = 1:numel(msh.f.ic{i,1})
                        k = msh.c.f.if(msh.f.ic{i,1}(j),:) == i;
                        for l = ["a","x"]
                            if ismembc(i,iu.(l))
                                for m = 1:size(s.c.x.x.xcunc.(l),3)
                                    s.c.x.x.xcunc.(l)(msh.f.ic{i,1}(j),k,m) = msh.c.f.Nf(msh.f.ic{i,1}(j),k)*s.c.x.x.xfunc.(l)(i,m)';
                                end
                            end
                        end
                    end
                end
                %  > dx (w/ analytical solution)...
                s.c.x.x.dx.a(:,:,1) = s.c.x.x.xcunc.a(:,:,1)-s.c.x.x.xcunc.x(:,:,1);            %  > a0+a1*h^n+...
                s.c.x.x.dx.a(:,:,2) = f.st                  -s.c.x.x.xcunc.x(:,:,1);            %  > a0.
                s.c.x.x.dx.a(:,:,3) = s.c.x.x.xcunc.a(:,:,1)-f.st;                              %  > a1*h^n+...
                %  > dx (w/ discrete solution)...
                if inp.PS.EE
                    if ~any(any(cellfun(@isempty,s.r.s.M),1),2)
                        %  > ...for tau-extrapolation.
                        s.c.x.x.dx.x{1}(:,:,1) = s.c.x.x.xcunc.x(:,:,1)-s.c.x.x.xcunc.x(:,:,2); %  > a1*h^n.
                        %  > ...for HO estimate.
                        s.c.x.x.dx.x{2}(:,:,1) = s.c.x.x.xcunc.x(:,:,3)-s.c.x.x.xcunc.x(:,:,1); %  > a0+a1*h^n.
                        s.c.x.x.dx.x{2}(:,:,2) = s.c.x.x.xcunc.x(:,:,4)-s.c.x.x.xcunc.x(:,:,1); %  > a0.
                        s.c.x.x.dx.x{2}(:,:,3) = s.c.x.x.xcunc.x(:,:,3)-s.c.x.x.xcunc.x(:,:,4); %  > a1*h^n.
                    end
                end
            end
        end
        %  > 1.5.1. -------------------------------------------------------
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
        %  > 1.5.2. -------------------------------------------------------
        %  > Auxiliary function #2.
        function [xfB] = xf_i(iu,s,vf)
            %  > tfunc*vf.
            xfA = cellfun(@(x,y) x*y,s.x.s.tfunc(iu,1),vf(iu,1),"UniformOutput",false);
            %  > Reshape...
            xfB = cat(2,xfA{:,1})';
        end
    end
end