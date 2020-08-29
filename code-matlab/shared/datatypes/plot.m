function h_out = plot(varargin)
% overload plot function to accept ts, tsd, iv structs
%
% MvdM 2016-01-20

if isstruct(varargin{1})
    if isfield(varargin{1},'type')
        switch varargin{1}.type
            case 'tsd'
                if size(varargin{1}.data,1) == 1
                    h = builtin('plot',varargin{1}.tvec,varargin{1}.data,varargin{2:end});
                elseif size(varargin{1}.data,1) == 2
                    h = builtin('plot',varargin{1}.data(1,:),varargin{1}.data(2,:),varargin{2:end});
                else
                   error('plotting for tsds with more than two dimensions not implemented yet.') 
                end
            case 'ts'
                cfg = []; cfg.openNewFig = 0;
                h = MultiRaster(cfg,varargin{1});
            otherwise
                 error('%s plot not implemented yet.',varargin{1}.type);   
        end
    else
        error('Unknown or undefined struct type.');
    end
else
    h = builtin('plot',varargin{:});
end

if nargout == 1
    h_out = h;
end