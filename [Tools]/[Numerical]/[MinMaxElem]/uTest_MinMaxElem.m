function uTest_MinMaxElem(doSpeed)
% Automatic test: MinMaxElem
% This is a routine for automatic testing. It is not needed for processing and
% can be deleted or moved to a folder, where it does not bother.
%
% uTest_MinMaxElem(doSpeed)
% INPUT:
%   doSpeed: If FALSE, the speed is not tested. Optional.
% OUTPUT:
%   On failure the test stops with an error.
%
% Tested: Matlab 6.5, 7.7, 7.8, WinXP
% Author: Jan Simon, Heidelberg, (C) 2009-2011 matlab.THISYEAR(a)nMINUSsimon.de

% $JRev: R-k V:037 Sum:3kml2816X1tm Date:04-Apr-2011 02:41:26 $
% $License: NOT_RELEASED $
% $File: Published\MinMaxElem\uTest_MinMaxElem.m $

% Initialize: ==================================================================
FuncName = 'uTest_MinMaxElem';  % $Managed by AutoFuncPath$

nRandTest = 5000;
nRandDisp = 1000;   % For showing dots as progress display

% Matlab 6.5 bug:
%   min(single([-Inf, NaN]))      =>  NaN instead of -Inf
%   min(single([-Inf, NaN, -5]))  =>  -5
bugMatlab65 = round(sscanf(version, '%f', 1) * 10) == 65;

disp(['== Test MinMaxElem:  ', datestr(now, 0)]);
disp(['Versions: ', char(10), ...
   which('MinMaxElem'), char(10)]);
pause(0.01);

if nargin == 0
   doSpeed = true;
end
minT = eps;

IntClassList = {'int8', 'uint8', 'int16', 'uint16', 'int32', 'uint32', ...
   'int64', 'uint64'};
ClassList    = cat(2, {'double', 'single'}, IntClassList);

for iClass = 1:length(ClassList)
   aClass  = ClassList{iClass};

   isFloat = strcmpi(aClass, 'single') || strcmpi(aClass, 'double');
   isInt64 = any(strfind(aClass, '64'));
   
   fprintf('== Known answer tests (%s)\n', upper(aClass));
   
   aZero = mycast(0, aClass);
   aOne  = mycast(1, aClass);
   aTwo  = mycast(2, aClass);
   
   q = mycast([], aClass);
   [x, y] = MinMaxElem(q);
   if isempty(x) && isempty(y)
      HideDisp('  ok: []');
   else
      error([FuncName, ': MinMaxElem([]): answers not empty']);
   end
   
   [x, y] = MinMaxElem(q, q);
   if isempty(x) && isempty(y)
      HideDisp('  ok: [], []');
   else
      error([FuncName, ': MinMaxElem([], []): answers not empty']);
   end
   
   [x, y] = MinMaxElem(q, q, q);
   if isempty(x) && isempty(y)
      HideDisp('  ok: [], [], []');
   else
      error([FuncName ': MinMaxElem([], [], []): answers not empty']);
   end
   
   [x, y] = MinMaxElem(aOne, q);
   if isequal(x, aOne) && isequal(y, aOne)
      HideDisp('  ok: 1, []');
   else
      error([FuncName, ': MinMaxElem(1, []): wrong answer: ', ...
         MatlabBugsprintfGG(x, y)]);
   end
   
   [x, y] = MinMaxElem(aOne);
   if isequal(x, aOne) && isequal(y, aOne)
      HideDisp('  ok: 1');
   else
      error([FuncName, ': MinMaxElem(1): wrong answer: ', ...
         MatlabBugsprintfGG(x, y)]);
   end
   
   [x, y] = MinMaxElem(aOne, aOne);
   if isequal(x, aOne) && isequal(y, aOne)
      HideDisp('  ok: 1, 1');
   else
      error([FuncName, ': MinMaxElem(1, 1): wrong answer: ', ...
         MatlabBugsprintfGG(x, y)]);
   end
   
   [x, y] = MinMaxElem(aZero, aOne);
   if isequal(x, aZero) && isequal(y, aOne)
      HideDisp('  ok: 0, 1');
   else
      error([FuncName, ': MinMaxElem(0, 1): wrong answer: ', ...
         MatlabBugsprintfGG(x, y)]);
   end
   
   [x, y] = MinMaxElem(aOne, aZero);
   if isequal(x, aZero) && isequal(y, aOne)
      HideDisp('  ok: 1, 0');
   else
      error([FuncName, ': MinMaxElem(1, 0): wrong answer: ', ...
         MatlabBugsprintfGG(x, y)]);
   end
   
   [x, y] = MinMaxElem([aOne, aZero]);
   if isequal(x, aZero) && isequal(y, aOne)
      HideDisp('  ok: [1, 0]');
   else
      error([FuncName, ': MinMaxElem([1, 0]): wrong answer: ', ...
         MatlabBugsprintfGG(x, y)]);
   end
   
   [x, y] = MinMaxElem([aZero, aZero]);
   if isequal(x, aZero) && isequal(y, aZero)
      HideDisp('  ok: [0, 0]');
   else
      error([FuncName, ': MinMaxElem([0, 0]): wrong answer: ', ...
         MatlabBugsprintfGG(x, y)]);
   end
   
   [x, y] = MinMaxElem([aOne, aZero], q);
   if isequal(x, aZero) && isequal(y, aOne)
      HideDisp('  ok: [1, 0]');
   else
      error([FuncName, ': MinMaxElem([1, 0], []): wrong answer: ', ...
         MatlabBugsprintfGG(x, y)]);
   end
   
   % NaN and Inf tests for SINGLE and DOUBLE only:
   if isFloat
      aInf = mycast(Inf, aClass);
      mInf = mycast(-Inf, aClass);
      aNaN = mycast(NaN, aClass);
      
      [x, y] = MinMaxElem([aInf, aZero]);
      if isequal(x, aZero) && isequal(y, aInf)
         HideDisp('  ok: [inf, 0]');
      else
         error([FuncName, ': MinMaxElem([inf, 0]): wrong answer: ', ...
            MatlabBugsprintfGG(x, y)]);
      end
      
      [x, y] = MinMaxElem([aOne, mInf]);
      if isequal(x, mInf) && isequal(y, aOne)
         HideDisp('  ok: [1, -inf]');
      else
         error([FuncName, ': MinMaxElem([1, -inf]): wrong answer: ', ...
            MatlabBugsprintfGG(x, y)]);
      end
      
      [x, y] = MinMaxElem(aNaN);
      if isempty(x) && isempty(y)
         HideDisp('  ok: NaN');
      else
         error([FuncName, ': MinMaxElem(NaN): wrong answer: ', ...
            MatlabBugsprintfGG(x, y)]);
      end
      
      [x, y] = MinMaxElem([aOne, aNaN]);
      if isequal(x, aOne) && isequal(y, aOne)
         HideDisp('  ok: [1, NaN]');
      else
         error([FuncName, ': MinMaxElem([1, NaN]): wrong answer: ', ...
            MatlabBugsprintfGG(x, y)]);
      end
      
      [x, y] = MinMaxElem([aNaN, aZero]);
      if isequal(x, aZero) && isequal(y, aZero)
         HideDisp('  ok: [NaN, 0]');
      else
         error([FuncName, ': MinMaxElem([NaN, 0]): wrong answer: ', ...
            MatlabBugsprintfGG(x, y)]);
      end
      
      [x, y] = MinMaxElem([aNaN, aOne, aNaN]);
      if isequal(x, aOne) && isequal(y, aOne)
         HideDisp('  ok: [NaN, 1, NaN]');
      else
         error([FuncName, ': MinMaxElem([NaN, 1, NaN]): wrong answer: ', ...
            MatlabBugsprintfGG(x, y)]);
      end
      
      [x, y] = MinMaxElem([aNaN, aOne, aNaN], [aNaN, mycast(10, aClass)]);
      if isequal(x, 1) && isequal(y, 10)
         HideDisp('  ok: [NaN, 1, NaN], [NaN, 10]');
      else
         error([FuncName, ': MinMaxElem([NaN, 1, NaN], [NaN, 10]): ', ...
            'wrong answer: ', MatlabBugsprintfGG(x, y)]);
      end
      
      [x, y] = MinMaxElem([aNaN, aNaN, aOne], aNaN, aInf);
      if isequal(x, aOne) && isequal(y, aInf)
         HideDisp('  ok: [NaN, NaN, 1], NaN, Inf');
      else
         error([FuncName, ': MinMaxElem([NaN, NaN, 1], NaN, Inf): ', ...
            'wrong answer: ', MatlabBugsprintfGG(x, y)]);
      end
      
      [x, y] = MinMaxElem([aNaN, aNaN, aOne], [aNaN, aInf, aTwo]);
      if isequal(x, aOne) && isequal(y, aInf)
         HideDisp('  ok: [NaN, NaN, 1], [NaN, Inf, 2]');
      else
         error([FuncName, ': MinMaxElem([NaN, NaN, 1], [NaN, Inf, 2]): ', ...
            'wrong answer: ', MatlabBugsprintfGG(x, y)]);
      end
      
      [x, y] = MinMaxElem([aNaN, aInf, aNaN, aOne], [aNaN, aInf, aTwo]);
      if isequal(x, 1) && isequal(y, Inf)
         HideDisp('  ok: [NaN, Inf, NaN, 1], [NaN, Inf, 2]');
      else
         error([FuncName, ...
            ': MinMaxElem([NaN, Inf, NaN, 1], [NaN, Inf, 2]): ', ...
            'wrong answer: ', MatlabBugsprintfGG(x, y)]);
      end
      
      [x, y] = MinMaxElem([aNaN, mInf, aNaN, aOne], [aNaN, aInf, aTwo]);
      if isequal(x, mInf) && isequal(y, aInf)
         HideDisp('  ok: [NaN, -Inf, NaN, 1], [NaN, Inf, 2]');
      else
         error([FuncName, ...
            ': MinMaxElem([NaN, -Inf, NaN, 1], [NaN, Inf, 2]): ', ...
            'wrong answer: ', MatlabBugsprintfGG(x, y)]);
      end
   end
   
   disp('  ok: known answer tests');
   drawnow;
   
   % Now a rude kind of searching bugs: Shotgun debugging
   % Although this is no reliable technique, it could reveal exceptions I
   % have not thought of.
   if isFloat
      fprintf('  %d tests, random numbers of [0:127], NaNs and Infs: ', ...
         nRandTest);
   else
      fprintf('  %d tests, random numbers of [0:127]: ', nRandTest);
   end
   
   for i = 1:nRandTest
      if mod(i, nRandDisp) == 0
         fprintf(1, '.');
      end
      q = mycast(floor(rand(1, 10) * 127), aClass);
      
      if isFloat
         q(rand(1, 10) < .05) = aNaN;
         q(rand(1, 10) < .05) = aInf;
         q(rand(1, 10) < .05) = mInf;
      end
      
      [x, y] = MinMaxElem(q);
      if bugMatlab65 && (strcmp(aClass, 'single') || isInt64)
         % Matlab 6.5 MIN(Single) is buggy, MIN(INT64) not implemented!
         xx = mycast(min(double(q)), aClass);
         yy = mycast(max(double(q)), aClass);
      else
         xx = min(q);
         yy = max(q);
      end
      if ~isequal(x, xx) || ~isequal(y, yy)
         fprintf('\n');
         disp(q);
         error([FuncName, ': MinMaxElem differs from MIN / MAX']);
      end
   end
   
   disp(' ok');
   drawnow;
   
   % Shotgun debugging:
   fprintf('  %d tests, several fields: ', nRandTest);
   for i = 1:nRandTest
      if mod(i, nRandDisp) == 0
         fprintf(1, '.');
      end
      
      qc = cell(1, 10);
      for j = 1:10
         n = floor(rand(1,1) * 10);
         if n
            if bugMatlab65 && isInt64
               % Matlab 6.5 bug: isequal(int64(-1), int64(-1)) == 0!
               q = mycast(floor(rand(1, n) * 10), aClass);
            else
               q = mycast(floor(rand(1, n) * 10) - 5, aClass);
            end
            
            if isFloat
               q(rand(1, n) < 0.1) = aNaN;
               q(rand(1, n) < 0.1) = aInf;
               q(rand(1, n) < 0.1) = mInf;
            end
         else
            q = mycast([], aClass);
         end
         qc{j} = q;
      end
      [x, y] = MinMaxElem(qc{:});
      
      qv = cat(2, qc{:});
      if bugMatlab65 && (strcmp(aClass, 'single') || isInt64)
         xx = mycast(min(double(qv(:))), aClass);
         yy = mycast(max(double(qv(:))), aClass);
      else
         xx = min(qv(:));
         yy = max(qv(:));
      end
      if ~isequal(x, xx) || ~isequal(y, yy)
         fprintf('\n');
         disp(q);
         error([FuncName, ': MinMaxElem differs from MIN / MAX']);
      end
   end
   disp(' ok');
   
   % ---------------------------------------------------------------------------
   % The 'finite' flag has no influence for integer classes!
   if isFloat
      fprintf('== Known answer tests (%s, finite):\n', upper(aClass));
      
      q      = mycast([], aClass);
      [x, y] = MinMaxElem(q, 'finite');
      if isempty(x) && isempty(y)
         HideDisp('  ok: []');
      else
         error([FuncName, ': MinMaxElem([], finite): answers not empty']);
      end
      
      [x, y] = MinMaxElem(q, q, 'finite');
      if isempty(x) && isempty(y)
         HideDisp('  ok: [], []');
      else
         error([FuncName, ': MinMaxElem([], [], finite): answers not empty']);
      end
      
      [x, y] = MinMaxElem(q, q, q, 'finite');
      if isempty(x) && isempty(y)
         HideDisp('  ok: [], [], []');
      else
         error([FuncName, ...
            ': MinMaxElem([], [], [], finite): answers not empty']);
      end
      
      [x, y] = MinMaxElem(q, aOne, q, 'finite');
      if isequal(x, aOne) && isequal(y, aOne)
         HideDisp('  ok: 1, []');
      else
         error([FuncName, ': MinMaxElem(1, [], finite): wrong answer: ', ...
            MatlabBugsprintfGG(x, y)]);
      end
      
      [x, y] = MinMaxElem(aOne, 'finite');
      if isequal(x, aOne) && isequal(y, aOne)
         HideDisp('  ok: 1');
      else
         error([FuncName, ': MinMaxElem(1, finite): wrong answer: ', ...
            MatlabBugsprintfGG(x, y)]);
      end
      
      [x, y] = MinMaxElem(aOne, aOne, 'finite');
      if isequal(x, 1) && isequal(y, aOne)
         HideDisp('  ok: 1, 1');
      else
         error([FuncName, ': MinMaxElem(1, 1, finite): wrong answer: ', ...
            MatlabBugsprintfGG(x, y)]);
      end
      
      [x, y] = MinMaxElem(aZero, aOne, 'finite');
      if isequal(x, aZero) && isequal(y, aOne)
         HideDisp('  ok: 0, 1');
      else
         error([FuncName, ': MinMaxElem(0, 1, finite): wrong answer: ', ...
            MatlabBugsprintfGG(x, y)]);
      end
      
      [x, y] = MinMaxElem(aOne, aZero, 'finite');
      if isequal(x, aZero) && isequal(y, aOne)
         HideDisp('  ok: 1, 0');
      else
         error([FuncName, ': MinMaxElem(1, 0, finite): wrong answer: ', ...
            MatlabBugsprintfGG(x, y)]);
      end
      
      [x, y] = MinMaxElem([aOne, aZero], 'finite');
      if isequal(x, aZero) && isequal(y, aOne)
         HideDisp('  ok: [1, 0]');
      else
         error([FuncName, ': MinMaxElem([1, 0], finite): wrong answer: ', ...
            MatlabBugsprintfGG(x, y)]);
      end
      
      [x, y] = MinMaxElem([aZero, aZero], 'finite');
      if isequal(x, aZero) && isequal(y, aZero)
         HideDisp('  ok: [0, 0]');
      else
         error([FuncName, ': MinMaxElem([0, 0], finite): wrong answer: ', ...
            MatlabBugsprintfGG(x, y)]);
      end
      
      [x, y] = MinMaxElem([inf, aZero], 'finite');
      if isequal(x, aZero) && isequal(y, aZero)
         HideDisp('  ok: [inf, 0]');
      else
         error([FuncName, ': MinMaxElem([inf, 0], finite): wrong answer: ', ...
            MatlabBugsprintfGG(x, y)]);
      end
      
      [x, y] = MinMaxElem([aOne, mInf, aTwo], 'finite');
      if isequal(x, aOne) && isequal(y, aTwo)
         HideDisp('  ok: [1, -inf, 2]');
      else
         error([FuncName, ...
            ': MinMaxElem([1, -inf, 2], finite): wrong answer: ', ...
            MatlabBugsprintfGG(x, y)]);
      end
      
      [x, y] = MinMaxElem([aOne, aZero], q, 'finite');
      if isequal(x, aZero) && isequal(y, aOne)
         HideDisp('  ok: [1, 0]');
      else
         error([FuncName, ...
            ': MinMaxElem([1, 0], [], finite): wrong answer: ', ...
            MatlabBugsprintfGG(x, y)]);
      end
      
      [x, y] = MinMaxElem(mycast([5, 6], aClass), ...
         mycast([4, 7], aClass), mycast([3, 8], aClass), 'finite');
      if isequal(x, mycast(3, aClass)) && isequal(y, mycast(8, aClass))
         HideDisp('  ok: [5, 6], [4, 7], [3, 8]');
      else
         error([FuncName, ': MinMaxElem([5, 6], [4, 7], [3, 8], finite): ', ...
            'wrong answer: ', MatlabBugsprintfGG(x, y)]);
      end
      
      [x, y] = MinMaxElem(mycast(5:6, aClass), [mInf, Inf], ...
         mycast([3, 8], aClass), 'finite');
      if isequal(x, 3) && isequal(y, 8)
         HideDisp('  ok: [5, 6], [-Inf, Inf], [3, 8]');
      else
         error([FuncName, ...
            ': MinMaxElem([5, 6], [-Inf, Inf], [3, 8], finite): ', ...
            'wrong answer: ', MatlabBugsprintfGG(x, y)]);
      end
      
      [x, y] = MinMaxElem(NaN, 'finite');
      if isempty(x) && isempty(y)
         HideDisp('  ok: NaN');
      else
         error([FuncName, ': MinMaxElem(NaN, finite): wrong answer: ', ...
            MatlabBugsprintfGG(x, y)]);
      end
      
      [x, y] = MinMaxElem([aOne, aNaN], 'finite');
      if isequal(x, aOne) && isequal(y, aOne)
         HideDisp('  ok: [1, NaN]');
      else
         error([FuncName, ': MinMaxElem([1, NaN], finite): wrong answer: ', ...
            MatlabBugsprintfGG(x, y)]);
      end
      
      [x, y] = MinMaxElem([aNaN, aZero], 'finite');
      if isequal(x, aZero) && isequal(y, aZero)
         HideDisp('  ok: [NaN, 0]');
      else
         error([FuncName, ...
            ': MinMaxElem([NaN, 0], finite): wrong answer: ', ...
            MatlabBugsprintfGG(x, y)]);
      end
      
      [x, y] = MinMaxElem([aNaN, aOne, aNaN], 'finite');
      if isequal(x, aOne) && isequal(y, aOne)
         HideDisp('  ok: [NaN, 1, NaN]');
      else
         error([FuncName, ...
            ': MinMaxElem([NaN, 1, NaN], finite): wrong answer: ', ...
            MatlabBugsprintfGG(x, y)]);
      end
      
      [x, y] = MinMaxElem([aNaN, aOne, aNaN], [aNaN, mycast(10, aClass)], ...
         'finite');
      if isequal(x, aOne) && isequal(y, mycast(10, aClass))
         HideDisp('  ok: [NaN, 1, NaN], [NaN, 10]');
      else
         error([FuncName, ...
            ': MinMaxElem([NaN, 1, NaN], [NaN, 10], finite): ', ...
            'wrong answer: ', MatlabBugsprintfGG(x, y)]);
      end
      
      [x, y] = MinMaxElem([aNaN, aNaN, aOne], mInf, ...
         mycast(6, aClass), mycast(5, aClass), aInf, 'finite');
      if isequal(x, aOne) && isequal(y, mycast(6, aClass))
         HideDisp('  ok: [NaN, NaN, 1], NaN, Inf');
      else
         error([FuncName, ...
            ': MinMaxElem([NaN, NaN, 1], NaN, Inf, finite): ', ...
            'wrong answer: ', MatlabBugsprintfGG(x, y)]);
      end
      
      [x, y] = MinMaxElem([aNaN, aNaN, aOne], aNaN, aInf, 'finite');
      if isequal(x, aOne) && isequal(y, aOne)
         HideDisp('  ok: [NaN, NaN, 1], NaN, Inf');
      else
         error([FuncName, ': MinMaxElem([NaN, NaN, 1], NaN, Inf, finite): ', ...
            'wrong answer: ', MatlabBugsprintfGG(x, y)]);
      end
      
      [x, y] = MinMaxElem([aNaN, mycast(5:6, aClass)], ...
         [mInf, aInf], mycast([3, 8], aClass), 'finite');
      if isequal(x, 3) && isequal(y, 8)
         HideDisp('  ok: [NaN, 5, 6], [-Inf, Inf], [3, 8]');
      else
         error([FuncName, ': MinMaxElem([NaN, 5, 6], [-Inf, Inf], [3, 8], ', ...
            'finite): wrong answer: ', MatlabBugsprintfGG(x, y)]);
      end
      
      [x, y] = MinMaxElem([aNaN, mycast(5:6, aClass)], ...
         [mInf, aNaN, mycast(7, aClass)], ...
         mycast([3, 8], aClass), 'finite');
      if isequal(x, 3) && isequal(y, 8)
         HideDisp('  ok: [NaN, 5, 6], [-Inf, NaN, 7], [3, 8]');
      else
         error([FuncName, ': MinMaxElem([NaN, 5, 6], [-Inf, NaN, 7], ', ...
            '[3, 8], finite): wrong answer: ', MatlabBugsprintfGG(x, y)]);
      end
      
      [x, y] = MinMaxElem(aNaN, [mInf, mycast(7, aClass)], ...
         mycast([3, 8], aClass), 'finite');
      if isequal(x, 3) && isequal(y, 8)
         HideDisp('  ok: NaN, [-Inf, 7], [3, 8]');
      else
         error([FuncName, ': MinMaxElem(NaN, [-Inf, 7], [3, 8], finite): ', ...
            'wrong answer: ', MatlabBugsprintfGG(x, y)]);
      end
      
      [x, y] = MinMaxElem(mycast(1:6, aClass), mInf, aInf, aNaN, 'finite');
      if isequal(x, 1) && isequal(y, 6)
         HideDisp('  ok: 1:6, -Inf, Inf, NaN');
      else
         error([FuncName, ': MinMaxElem(1:6, -Inf, Inf, NaN, finite): ', ...
            'wrong answer: ', MatlabBugsprintfGG(x, y)]);
      end
      disp('  ok: known answer tests');
      drawnow;
      
      % Shotgun debugging:
      fprintf('  %d tests, random numbers, NaNs and Infs: ', ...
         nRandTest);
      for i = 1:nRandTest
         if mod(i, nRandDisp) == 0
            fprintf(1, '.');
         end
         q = mycast(floor(rand(1, 10) * 10) - 5, aClass);
         
         % The 'finite' flag is checked for floating point types only!
         % if isFloat
         q(rand(1, 10) < 0.05) = aNaN;
         q(rand(1, 10) < 0.05) = aInf;
         q(rand(1, 10) < 0.05) = mInf;
         % end
         
         [x, y]   = MinMaxElem(q, 'finite');
         [xx, yy] = localMinMaxFin(q);
         if ~isequal(x, double(xx)) || ~isequal(y, double(yy))
            fprintf('\n');
            disp(q);
            error([FuncName, ...
               ': MinMaxElem(finite) differs from local Matlab function.']);
         end
      end
      disp(' ok');
      drawnow;
      
      % Shotgun debugging:
      fprintf('  %d tests, several fields: ', nRandTest);
      for i = 1:nRandTest
         if mod(i, nRandDisp) == 0
            fprintf(1, '.');
         end
         
         qc = cell(1, 10);
         for j = 1:10
            n = floor(rand(1, 1) * 10);
            if n
               q = mycast(floor(rand(1, n) * 10), aClass);
               % The 'finite' flag is checked for floating point types only!
               % if isFloat
               q(rand(1, n) < 0.1) = aNaN;
               q(rand(1, n) < 0.1) = aInf;
               q(rand(1, n) < 0.1) = mInf;
               % end
            else
               q = mycast([], aClass);
            end
            qc{j} = q;
         end
         [x, y]   = MinMaxElem(qc{:}, 'finite');
         [xx, yy] = localMinMaxFin(qc{:});
         if ~isequal(x, double(xx)) || ~isequal(y, double(yy))
            fprintf('\n');
            disp(q);
            error([FuncName, ...
               ': MinMaxElem(finite) differs from local Matlab function.']);
         end
      end
      fprintf(' ok\n');
   end
   
   clear('aInf', 'mInf', 'aNaN');
end

% Special tests: ---------------------------------------------------------------
disp('== Test indices and argument number:');

[x, y, a, b, c, d] = MinMaxElem([]);
if isequal({x, y, a, b, c, d}, {[], [], [], [], [], []})
   HideDisp('  ok: 1');
else
   error([FuncName, ': MinMaxElem(1): wrong answer']);
end

[x, y, a, b, c, d] = MinMaxElem(1);
if isequal({x, y, a, b, c, d}, {1, 1, 1, 1, 1, 1})
   HideDisp('  ok: 1');
else
   error([FuncName, ': MinMaxElem(1): wrong answer']);
end

[x, y, a, b, c, d] = MinMaxElem(1:2);
if isequal({x, y, a, b, c, d}, {1, 2, 1, 2, 1, 1})
   HideDisp('  ok: 1:2');
else
   error([FuncName, ': MinMaxElem(1:2): wrong answer']);
end

[x, y, a, b, c, d] = MinMaxElem([2, 1]);
if isequal({x, y, a, b, c, d}, {1, 2, 2, 1, 1, 1})
   HideDisp('  ok: [2,1]');
else
   error([FuncName, ': MinMaxElem([2, 1]): wrong answer']);
end

[x, y, a, b, c, d] = MinMaxElem(1, 2);
if isequal({x, y, a, b, c, d}, {1, 2, 1, 1, 1, 2})
   HideDisp('  ok: 1, 2');
else
   error([FuncName, ': MinMaxElem(1, 2): wrong answer']);
end

[x, y, a, b, c, d] = MinMaxElem(3, 2, 1);
if isequal({x, y, a, b, c, d}, {1, 3, 1, 1, 3, 1})
   HideDisp('  ok: 3, 2, 1');
else
   error([FuncName, ': MinMaxElem(3, 2, 1): wrong answer']);
end


[x, y, a, b, c, d] = MinMaxElem(3:4, [2,1], [6,5]);
if isequal({x, y, a, b, c, d}, {1, 6, 2, 1, 2, 3})
   HideDisp('  ok: 3:4, [2,1], [6,5]');
else
   error([FuncName, ': MinMaxElem(3:4, [2,1], [6,5]): wrong answer']);
end

[x, y, a, b, c, d] = MinMaxElem(1:3, 2:4, 1:5);
if isequal({x, y, a, b, c, d}, {1, 5, 1, 5, 1, 3})
   HideDisp('  ok: 1:3, 2:4, 1:5');
else
   error([FuncName, ': MinMaxElem(1:3, 2:4, 1:5): wrong answer']);
end

[x, y, a, b, c, d] = MinMaxElem(2:6, 2:4, 1:5);
if isequal({x, y, a, b, c, d}, {1, 6, 1, 5, 3, 1})
   HideDisp('  ok: 2:6, 2:4, 1:5');
else
   error([FuncName, ': MinMaxElem(2:6, 2:4, 1:5): wrong answer']);
end

disp('  ok');

% Speed: -----------------------------------------------------------------------
% This test works for INT types also: MinMaxElem is 10-20% faster for small
% arrays, 50% faster for [1 x 1e3] and larger.
for iClass = 1:2
   aClass = ClassList{iClass};
   fprintf('== Test Speed (%s):\n', upper(aClass));
   
   % Find a suiting number of loops:
   q = mycast(1:1000, aClass);
   if doSpeed
      iLoop     = 0;
      startTime = cputime;
      while cputime - startTime < 1.0
         a = min(q); %#ok<*NASGU>
         clear('a');
         iLoop = iLoop + 1;
      end
      nLoops = 1000 * ceil(iLoop / ((cputime - startTime) * 1000));
      disp([sprintf('  %d', nLoops) ' loops on this machine.']);
   else
      disp('  Use at least 2 loops (displayed times are strange!)');
      nLoops = 2;
   end
   disp('  (*): Extrapolated time for slow MIN&MAX');
   
   q = mycast(1:100, aClass);
   tic;
   for i = 1:nLoops
      x = min(q);
      y = max(q);
      clear('x', 'y');
   end
   et = minT + toc;
   
   tic;
   for i = 1:nLoops
      [x, y] = MinMaxElem(q);
      clear('x', 'y');
   end
   et2 = toc;
   disp(['  [1:100]:        MIN & MAX:', sprintf('%6.2f', et), ...
      '   MinMaxElem:', sprintf('%6.2f', et2), '   ==> ', ...
      sprintf('%5.1f', 100 * et2 / et), '%']);
   
   q = mycast(1:1000, aClass);
   tic;
   for i = 1:nLoops
      x = min(q);
      y = max(q);
      clear('x', 'y');
   end
   et = minT + toc;
   
   tic;
   for i = 1:nLoops
      [x, y] = MinMaxElem(q);
      clear('x', 'y');
   end
   et2 = toc;
   disp(['  [1:1000]:       MIN & MAX:', sprintf('%6.2f', et), ...
      '   MinMaxElem:', sprintf('%6.2f', et2), '   ==> ', ...
      sprintf('%5.1f', 100 * et2 / et), '%']);
   
   q = mycast(1:10000, aClass);
   mLoops = ceil(nLoops / 10);
   tic;
   for i = 1:mLoops
      x = min(q);
      y = max(q);
      clear('x', 'y');
   end
   et = minT + nLoops * toc / mLoops;  % Extrapolate time
   
   tic;
   for i = 1:nLoops
      [x, y] = MinMaxElem(q);
      clear('x', 'y');
   end
   et2 = toc;
   disp(['  [1:10000]:      MIN & MAX:', sprintf('%6.2f', et), ...
      '   MinMaxElem:', sprintf('%6.2f', et2), '   ==> ', ...
      sprintf('%5.1f', 100 * et2 / et), '%  (*)']);
   
   q = mycast(1000:-1:1, aClass);
   tic;
   for i = 1:nLoops
      x = min(q);
      y = max(q);
      clear('x', 'y');
   end
   et = minT + toc;
   
   tic;
   for i = 1:nLoops
      [x, y] = MinMaxElem(q);
      clear('x', 'y');
   end
   et2 = toc;
   disp(['  [1000:-1:1]:    MIN & MAX:', sprintf('%6.2f', et), ...
      '   MinMaxElem:', sprintf('%6.2f', et2), '   ==> ', ...
      sprintf('%5.1f', 100 * et2 / et), '%']);
   
   q = mycast(rand(1, 100), aClass);
   tic;
   for i = 1:nLoops
      x = min(q);
      y = max(q);
      clear('x', 'y');
   end
   et = minT + toc;
   
   tic;
   for i = 1:nLoops
      [x, y] = MinMaxElem(q);
      clear('x', 'y');
   end
   et2 = toc;
   disp(['  rand(1, 100):   MIN & MAX:', sprintf('%6.2f', et), ...
      '   MinMaxElem:', sprintf('%6.2f', et2), '   ==> ', ...
      sprintf('%5.1f', 100 * et2 / et), '%']);
   
   q = mycast(rand(1, 1000), aClass);
   tic;
   for i = 1:nLoops
      x = min(q);
      y = max(q);
      clear('x', 'y');
   end
   et = minT + toc;
   
   tic;
   for i = 1:nLoops
      [x, y] = MinMaxElem(q);
      clear('x', 'y');
   end
   et2 = toc;
   disp(['  rand(1, 1000):  MIN & MAX:', sprintf('%6.2f', et), ...
      '   MinMaxElem:', sprintf('%6.2f', et2), '   ==> ', ...
      sprintf('%5.1f', 100 * et2 / et), '%']);
   
   q = mycast(rand(1, 10000), aClass);
   mLoops = ceil(nLoops / 10);
   tic;
   for i = 1:mLoops
      x = min(q);
      y = max(q);
      clear('x', 'y');
   end
   et = minT + nLoops * toc / mLoops;
   
   tic;
   for i = 1:nLoops
      [x, y] = MinMaxElem(q);
      clear('x', 'y');
   end
   et2 = toc;
   disp(['  rand(1, 10000): MIN & MAX:', sprintf('%6.2f', et), ...
      '   MinMaxElem:', sprintf('%6.2f', et2), '   ==> ', ...
      sprintf('%5.1f', 100 * et2 / et), '%  (*)']);
   
   q1 = mycast(rand(100), aClass);
   q2 = mycast(rand(100), aClass);
   tic;
   for i = 1:nLoops
      x = min(min(q1(:)), min(q2(:)));
      y = max(max(q1(:)), max(q2(:)));
      clear('x', 'y');
   end
   et = minT + toc;
   
   tic;
   for i = 1:nLoops
      [x, y] = MinMaxElem(q1, q2);
      clear('x', 'y');
   end
   et2 = toc;
   disp(['  2 x rand(100):  MIN & MAX:', sprintf('%6.2f', et), ...
      '   MinMaxElem:', sprintf('%6.2f', et2), '   ==> ', ...
      sprintf('%5.1f', 100 * et2 / et), '%']);
   
   q = mycast(rand(1, 1E7), aClass);
   mLoops = ceil(nLoops / 3000);
   tic;
   for i = 1:mLoops
      x = min(q);
      y = max(q);
      clear('x', 'y');
   end
   et = minT + toc;
   
   tic;
   for i = 1:mLoops
      [x, y] = MinMaxElem(q);
      clear('x', 'y');
   end
   et2 = toc;
   disp(['  rand(1, 1E7):   MIN & MAX:', sprintf('%6.2f', et), ...
      '   MinMaxElem:', sprintf('%6.2f', et2), '   ==> ', ...
      sprintf('%5.1f%%  (%d loops)', 100 * et2 / et, mLoops)]);
   
   fprintf('\n');
end

fprintf('== MinMaxElem passed the tests\n');

return;

% ******************************************************************************
function Str = MatlabBugsprintfGG(x, y)
% sprintf('%g %g', [], 2) returns empty string in Matlab 6.5!

if isempty(x)
   Str = '[], ';
else
   Str = sprintf('%g ', double(x));
end

if isempty(y)
   Str = [Str, '[]'];
else
   Str = [Str, sprintf('%g ', double(y))];
end

return;

% ******************************************************************************
function A = mycast(A, ClassName)
% Simulate CAST for Matlab 6.5
A = feval(ClassName, A);
return;

% ******************************************************************************
function HideDisp(S)  %#ok<INUSD>
% Uncomment this to display the large number of successful checks
% disp(S);
return;

% ******************************************************************************
function [XMin, XMax] = localMinMaxFin(varargin)
% Find min and max finite element over one or several arrays

% Was JRev: R0b Build004 Sum:216F1F62065 Date:11-Dec-2003 17:48:55
% Was File: Tools\GLMath\MinMaxFin.m

% Do the work: =================================================================
XMin = [];
XMax = [];
for i = 1:nargin
   X = varargin{i}(:);
   X = X(isfinite(X));
   if length(X)
      if isempty(XMin)  % This must be considered
         XMin = min(X);
      else
         XMin = min(XMin, min(X));
      end
      
      if isempty(XMax)
         XMax = max(X);
      else
         XMax = max(XMax, max(X));
      end
   end
end

return;
