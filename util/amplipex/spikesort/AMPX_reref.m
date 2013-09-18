function data = AMPX_reref(data,av)
% function data_out = AMPX_reref(data_in,av)
%
% subtracts av from all channels in data_in
%
% NOTE: should add provenance tracking to this

for iCh = length(data.channels)
   
    if isa(data.channels{iCh},'int16')
       
        data.channels{iCh} = double(data.channels{iCh});
        
    end
    
    data.channels{iCh} = data.channels{iCh} - av;
    
end

data.rereferenced = 'yes';