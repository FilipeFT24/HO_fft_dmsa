classdef AuxiliaryFunctions
    methods (Static)            
        %% > 1. -----------------------------------------------------------
        % >> 1.1. ---------------------------------------------------------
        %  > Similar (faster) version of other built-in functions (to compute distance between 2 points).
        function [D] = dist(XY)
            D = sqrt((XY(1,1)-XY(2,1)).^2+(XY(1,2)-XY(2,2)).^2);
        end
        % >> 1.2. ---------------------------------------------------------
        %  > Similar (faster) version of built-in function "mean".
        function [y] = mean(x,j)
            y = sum(x,j,"default","includenan")./size(x,j);
        end
        % >> 1.3. ---------------------------------------------------------
        %  > Similar (faster) version of built-in function "setdiff" (w/o verification(s)).
        function [z] = setdiff(x,y)
            w    = false(1,max(max(x),max(y)));
            w(x) = true;
            w(y) = false;
            z    = x(w(x));
        end
        % >> 1.4. ---------------------------------------------------------
        %  > Similar (faster) version of following expression: arrayfun(@(x) find(a == x),b).
        function [j,k] = find_c(a,b)
            i = repmat(reshape(a,[],1),1,numel(b)) == reshape(b,1,[]); [j,~] = find(i);
            k = any(i,2);
        end
        % >> 1.5. ---------------------------------------------------------
        %  > Similar (faster) version of unique(v,'rows').
        function [v] = unique_r(u)
            [us,is] = sortrows(u);
            k       = find([true;any(diff(us,1,1),2);true]);
            js      = is(k(diff(k) >= 1));
            v       = u(js,:);
        end
        % >> 1.6. ---------------------------------------------------------
        %  > Compute error slope.
        function [s] = Slope(h,e)
            [m,n]  = size(e);
            i      = 1:n;
            j      = 1:m-1;
            s(j,i) = log(e(j+1,i)./e(j,i))./log(h(j+1)./h(j));
        end
    end
end