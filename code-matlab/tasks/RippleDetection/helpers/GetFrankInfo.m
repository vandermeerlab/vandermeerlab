function info_out = GetFrankInfo(info_in,fieldname,day,epoch,tetrode,varargin)
%GETFRANKINFO Retrive data contained in the terminal structs of Frank lab 
% cellinfo or tetinfo arrays 
%   info_out = GetFrankInfo(info_in,fieldname,day,epoch,tetrode)
%   info_out = GetFrankInfo(info_in,fieldname,day,epoch,tetrode,cell)
%
%   INPUTS
%       info_in:   Cellinfo or tetinfo cell array
%       fieldname: String specifying the data you want to retrieve, e.g.
%                  'area' or 'numspikes' etc.
%
%       The remainder are indices into the cell array:
%       day:        the recording day, e.g. 3
%       epoch:      The epoch, e.g. 7
%       tetrode:    The tetrode number, e.g. 29
%       varargin:   The cell (unit) number, e.g. 3
%
% aacarey Dec 2017

assignIfEmpty = NaN;

switch nargin
    case 5
        if ~isempty(info_in{1,day}) && ~isempty(info_in{1,day}{1,epoch}) && ~isempty(info_in{1,day}{1,epoch}{1,tetrode}) && isfield(info_in{1,day}{1,epoch}{1,tetrode},fieldname)
           info_out = info_in{1,day}{1,epoch}{1,tetrode}.(fieldname);
        else
           info_out = assignIfEmpty;
        end
    case 6
        if ~isempty(info_in{1,day}) && ~isempty(info_in{1,day}{1,epoch}) && ~isempty(info_in{1,day}{1,epoch}{1,tetrode}) && ~isempty(info_in{1,day}{1,epoch}{1,tetrode}{1,varargin{1}}) && isfield(info_in{1,day}{1,epoch}{1,tetrode}{1,varargin{1}},fieldname)
            info_out = info_in{1,day}{1,epoch}{1,tetrode}{1,varargin{1}}.(fieldname);
        else
            info_out = assignIfEmpty;
        end
    otherwise
        error('Unsupported number of input arguments')
end

end

