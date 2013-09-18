function raise_me (me,nothing)
% UltraMegaSort2000 by Hill DN, Mehta SB, & Kleinfeld D  - 07/12/2010
%
% raise_me - callback to raise clicked object to front of visibility
%
% Usage:
%      raise_me (hObject, event, me)
%
% Description:  
% Pulls the clicked object to the top of the ui stack; useful
% for raising partially masked objects to the front of a plot.
% GUI-shortcut for 'uistack':  Left-clicking brings to the top,
% right-clicking sends to bottom.
%
% Input: 
%   me      - handle to the current object
%   nothing - unused callback argument
%
[o,h] = gcbo;

switch(get(h, 'SelectionType')),
    case ('normal'),
        uistack(me, 'top');
    case ('alt'),
        uistack(me, 'bottom');
end

drawnow;
