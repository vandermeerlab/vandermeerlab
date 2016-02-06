function h_out = plot(varargin)
% overload plot function to accept ts, tsd, iv structs
%
% MvdM 2016-01-20

if isstruct(varargin{1})
    if isfield(varargin{1},'type')
        switch varargin{1}.type
            case 'tsd'
                h = builtin('plot',varargin{1}.tvec,varargin{1}.data,varargin{2:end});
            case 'ts'
                h = MultiRaster(varargin{1});
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