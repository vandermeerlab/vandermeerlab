function pass_flag = isFrankRat(ratID)
%ISFRANKRAT Determine if rat is a Frank rat
%   pass_flag = isFrankRat(ratID)
%
%   INPUTS
%       ratID - string specifying rat ID (ex: 'bon' or 'R050')
%
%   OUTPUTS
%       pass_flag - 1 if rat is a Frank rat, 0 if rat is not a Frank rat
%
%   NOTE: this function contains a list of approved Frank Rats. To add a 
%         new rat, open isFrankRat and edit the rat list. Ensure
%         that all data and folders are structured as expected by the
%         analysis software for the ripple detection project.
%
% aacarey Dec 2017

FrankRats = {'bon','cor'};

if contains(ratID,FrankRats)
    pass_flag = 1;
else
    pass_flag = 0;
end

end

