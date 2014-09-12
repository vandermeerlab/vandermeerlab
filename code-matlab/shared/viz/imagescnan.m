function h = imagescnan(varargin)
% IMAGESC with NaNs transparent
%
% MvdM 2014

switch (nargin)
  case 1
    hh = imagesc(varargin{1});
    a = varargin{1};
  case 3
    hh = imagesc(varargin{:});
    a = varargin{3};
end

set(hh,'alphadata',~isnan(a));

if nargout > 0
    h = hh;
end