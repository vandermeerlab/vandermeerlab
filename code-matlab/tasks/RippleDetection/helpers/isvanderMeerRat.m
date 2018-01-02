function pass_flag = isvanderMeerRat(ratID)
%ISVANDERMEERRAT Determine if rat is a van der Meer rat
%   pass_flag = isvanderMeerRat(ratID)
%
%   INPUTS
%       ratID - string specifying rat ID (ex: 'bon' or 'R050')
%
%   OUTPUTS
%       pass_flag - 1 if rat is a van der Meer rat, 0 if rat is not a van
%                   der Meer rat
%
%   NOTE: this function contains a list of approved van der Meer rats. To 
%         add a new rat, open this function and edit the rat list. Ensure
%         that all data and folders are structured as expected by the
%         analysis software for the ripple detection project.
%
% aacarey Dec 2017

vdMRats = {'R050','R064'};

if contains(ratID,vdMRats)
    pass_flag = 1;
else
    pass_flag = 0;
end

end

