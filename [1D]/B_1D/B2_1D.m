classdef B2_1D
    methods(Static)
        %% > 1. -----------------------------------------------------------
        % >> Set up 'P-Standard' and 'P-Adaptive' runs.
        function [obj] = Run_P(inp,el,t)
            %  > 'msh'.
            msh = A2_1D.Set_msh(el,t);
            switch inp.T
                case 1
                    % >> 'P-Standard' run.
                    %  > 'obj'.
                    obj.f = A3_1D.Initialize_f(inp,msh);
                    obj.s = B1_1D.Initialize_s(inp,msh);
                    obj.s = B1_1D.Update_s    (inp,msh,obj.f,obj.s);
                    obj.e = NE_XD.Initialize_e(inp,msh);
                    obj.e = NE_XD.Update_e    (    msh,obj.e,obj.s);
                    %  > Plot...
                    for i = reshape(setdiff(string(fieldnames(obj.s)),"u"),1,[])
                        if any(any(cellfun(@isempty,obj.s.(i).s.M),1),2)
                            return;
                        end
                    end
                    PP_1D.Plot(inp,msh,obj,...
                        {[1,0,0,0,0],...   %  > 1) 1D.
                         [0,0],...         %  > 2) AMR (AD #1).
                         [0],...           %  > 3) AMR (AD #2).      
                         [0,0],...         %  > 4) AMR (AD #3).
                         [0,0],...         %  > 5) AMR (NE).
                         [0],...           %  > 6) Ref. level.
                         [0],...           %  > 7) Appendix B.
                         [0,0,0,0]},{[]}); %  > 8) 1D.
                    %  > Print to terminal...
                    if ~inp.PS.EE
                        fprintf("  E(c)    = %.6e %.6e\n",obj(end).e.a.n_abs.c  (:,1));
                        fprintf("Tau(f) #k = %.6e %.6e\n",obj(end).e.a.n_abs.t.f(:,3));
                        fprintf("eta       = %.6f     \n",obj(end).e.a.n_abs.c(1,1)./obj(end).e.a.n_abs.c(2,1));
                        fprintf("NNZ(A)    = %d       \n",obj(end).s.c.m.nnz.At);
                        fprintf("h         = %.6e     \n",msh.h);
                        fprintf("NC        = %d       \n",msh.c.Nc);
                    else
                        fprintf("  E(c)    (A) = %.6e %.6e\n",obj(end).e.a   .n_abs.c  (:,1));
                        fprintf("          (l) = %.6e %.6e\n",obj(end).e.x(1).n_abs.c  (:,1));
                        fprintf("          (h) = %.6e %.6e\n",obj(end).e.x(2).n_abs.c  (:,1));
                        fprintf("Tau(f) #k (A) = %.6e %.6e\n",obj(end).e.a   .n_abs.t.f(:,3));
                        fprintf("          (l) = %.6e %.6e\n",obj(end).e.x(1).n_abs.t.f(:,1));
                        fprintf("          (h) = %.6e %.6e\n",obj(end).e.x(2).n_abs.t.f(:,3));
                        fprintf("eta           = %.6f     \n",obj(end).e.a.n_abs.c(1,1)./obj(end).e.a.n_abs.c(2,1));
                        fprintf("NNZ(A)        = %d       \n",obj(end).s.c.m.nnz.At);
                        fprintf("h             = %.6e     \n",msh.h);
                        fprintf("NC            = %d       \n",msh.c.Nc);
                    end
                case 2
                    % >> 'P-Adaptive' run.
                    obj = P_Ad.P_Adaptive(inp,msh);
                otherwise
                    return;
            end
        end
    end
end