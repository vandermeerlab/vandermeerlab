function varargout = ginput_ax(ha,n)
if nargin<2
    n=1;
end
k = 0;
button = 0;
xy = zeros(n,2);
hf = get(ha,'parent');
figure(hf);
set(hf,'WindowButtonMotionFcn',@changepointer)
set(ha,'ButtonDownFcn',@getpoints)
hp = get(ha,'children');
ht = get(hp,'hittest');
set(hp,'hittest','off')
axlim = get(ha,'Position');
fglim = get(hf,'Position');
x1 = axlim(1)*fglim(3) + fglim(1);
x2 = (axlim(1)+axlim(3))*fglim(3) + fglim(1);
y1 = axlim(2)*fglim(4) + fglim(2);
y2 = (axlim(2)+axlim(4))*fglim(4) + fglim(2);
waitfor(hf,'WindowButtonMotionFcn',[])
if iscell(ht)
    for jj=1:length(ht)
        set(hp(jj),'hittest',ht{jj})
    end
else
    set(hp,'hittest',ht)
end

% Mouse-Button recognition...
if(strcmp(button, 'normal'))
    button = 1; % left
elseif(strcmp(button, 'extend'))
    button = 2; % right
elseif(strcmp(button, 'alt'))
    button = 3; % middle
else
    button = 4; % double click any mousebutton
end

if nargout==3 
    varargout{1} = xy(:,1);
    varargout{2} = xy(:,2);
    varargout{3} = button;
elseif nargout==2
    varargout{1} = xy(:,1);
    varargout{2} = xy(:,2);
else
    varargout{1} = xy;
end
    function changepointer(~,~)
        pntr = get(0,'PointerLocation');
        if pntr(1)>x1 && pntr(1)<x2 && pntr(2)>y1 && pntr(2)<y2
            set(hf,'Pointer','crosshair')
        else
            set(hf,'Pointer','arrow')
        end
    end
    function getpoints(src,evnt)
        cp = get(src,'CurrentPoint');
        button = get(hf, 'SelectionType');
        k = k+1;
        xy(k,:) = cp(1,1:2);
        if k==n
            set(hf,'Pointer','arrow')
            set(hf,'WindowButtonMotionFcn',[])
            set(ha,'ButtonDownFcn',[])
        end
    end
end