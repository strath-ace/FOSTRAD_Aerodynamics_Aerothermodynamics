function isbin =  isbinary(theFile,nChars)
%ISBINARY   True for binary files
%
%   Syntax:
%      IS = ISBINARY(FILE,N)
%
%   Inputs:
%      FILE
%      N     Number of elements in file where char 0 will be searched
%            [ 1000 ]
%
%   Output:
%      IS   1 if FILE is binary and 0 otherwise, or [] if file fid=-1
%
%   MMA 17-09-2005, martinho@fis.ua.pt

%   Department of Physics
%   University of Aveiro, Portugal

%% License
% Copyright (c) 2009, M MA
% All rights reserved.
% 
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are met:
% 
% * Redistributions of source code must retain the above copyright notice, this
%   list of conditions and the following disclaimer.
% 
% * Redistributions in binary form must reproduce the above copyright notice,
%   this list of conditions and the following disclaimer in the documentation
%   and/or other materials provided with the distribution
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
% DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE
% FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
% DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
% SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
% CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
% OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
% OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
%%

if nargin < 2
  nChars = 1000;
end

isbin=0;
fid=fopen(theFile);
if fid~=-1
  c=fread(fid,nChars);
  if any(c==0), isbin = 1; end
  fclose(fid);
else
  isbin=[];
end
