function rat_list = SWRrats
%SWRRATS Get list of rat IDs for the SWR detection project
%   rat_list = SWRrats
%            = {'R050','R064','RBon','RCon'};
%
%   NOTE: this function contains a list of approved ripple project rats. To 
%         add a new rat, open this function and edit the rat list. Ensure
%         that all data and folders are structured as expected by the
%         analysis software for the ripple detection project.
%
% aacarey Dec 2017

mfun = mfilename;

rat_list = {'R050','R064','bon','cor'};
disp([mfun,': getting rat IDs for Digit, Obex, Bond, and Coriander...'])

end

