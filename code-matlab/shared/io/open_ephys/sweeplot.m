function varargout = sweeplot(varargin)
% MvdM evoked potential viewer for open ephys data
%
% add folder that contains this (.m and .fig) to path!
%
% user configurables:
%
% params.decimateFactor = 30;
% params.timeWindow = [-0.05 0.1]; % time window for viewer, centered on event times
% params.sweepChannel = 0; % I/O ID for sweep start signal
% params.eventChannel = 1; % I/O ID for stim events

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @sweeplot_OpeningFcn, ...
                   'gui_OutputFcn',  @sweeplot_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% SET PARAMETERS -- USER TO CONFIGURE %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global params;
params.decimateFactor = 30;
params.timeWindow = [-0.05 0.1]; % time window for viewer, centered on event times
params.sweepChannel = 0; % I/O ID for sweep start signal
params.eventChannel = 1; % I/O ID for stim events

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% END PARAMETERS -- DO NOT EDIT BELOW UNLESS KNOW WHAT DOING %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% --- Executes just before sweeplot is made visible.
function sweeplot_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to sweeplot (see VARARGIN)

% Choose default command line output for sweeplot
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);


% UIWAIT makes sweeplot wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = sweeplot_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

% --- Executes on button press in ReloadButton.
function ReloadButton_Callback(hObject, eventdata, handles)
% hObject    handle to ReloadButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
LoadData(handles)

% --------------------------------------------------------------------
function FileMenu_Callback(hObject, eventdata, handles)
% hObject    handle to FileMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function OpenMenuItem_Callback(hObject, eventdata, handles)
% hObject    handle to OpenMenuItem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
file = uigetfile('*.fig');
if ~isequal(file, 0)
    open(file);
end

% --------------------------------------------------------------------
function PrintMenuItem_Callback(hObject, eventdata, handles)
% hObject    handle to PrintMenuItem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
printdlg(handles.figure1)

% --------------------------------------------------------------------
function CloseMenuItem_Callback(hObject, eventdata, handles)
% hObject    handle to CloseMenuItem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
selection = questdlg(['Close ' get(handles.figure1,'Name') '?'],...
                     ['Close ' get(handles.figure1,'Name') '...'],...
                     'Yes','No','Yes');
if strcmp(selection,'No')
    return;
end

delete(handles.figure1)

function EventBox_Callback(hObject, eventdata, handles)
% hObject    handle to EventBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of EventBox as text
%        str2double(get(hObject,'String')) returns contents of EventBox as a double


% --- Executes during object creation, after setting all properties.
function EventBox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to EventBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in EventLeft.
function EventLeft_Callback(hObject, eventdata, handles)
% hObject    handle to EventLeft (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in EventRight.
function EventRight_Callback(hObject, eventdata, handles)
% hObject    handle to EventRight (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function SweepBox_Callback(hObject, eventdata, handles)
% hObject    handle to SweepBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of SweepBox as text
%        str2double(get(hObject,'String')) returns contents of SweepBox as a double

% find events for this sweep
FindEvents(handles);

% plot
UpdatePlot(handles);


% --- Executes during object creation, after setting all properties.
function SweepBox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to SweepBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in SweepLeft.
function SweepLeft_Callback(hObject, eventdata, handles)
% hObject    handle to SweepLeft (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
currentSweep = str2double(get(handles.SweepBox,'String'));
if currentSweep ~= 1 % last sweep
    set(handles.SweepBox,'String',num2str(currentSweep-1));
end

% find events for this sweep
FindEvents(handles);

% plot
UpdatePlot(handles);

% --- Executes on button press in SweepRight.
function SweepRight_Callback(hObject, eventdata, handles)
% hObject    handle to SweepRight (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global data;

currentSweep = str2double(get(handles.SweepBox,'String'));
if currentSweep ~= data.nSweeps; % last sweep
    set(handles.SweepBox,'String',num2str(currentSweep+1));
end

% find events for this sweep
FindEvents(handles);

% plot
UpdatePlot(handles);


function FileNameBox_Callback(hObject, eventdata, handles)
% hObject    handle to FileNameBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of FileNameBox as text
%        str2double(get(hObject,'String')) returns contents of FileNameBox as a double


% --- Executes during object creation, after setting all properties.
function FileNameBox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to FileNameBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in LoadButton.
function LoadButton_Callback(hObject, eventdata, handles)
% hObject    handle to LoadButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[input_file,pathname] = uigetfile( ...
       {'*.continuous', 'OpenEphys Data File (*.continuous)'; ...
        '*.*', 'All Files (*.*)'}, ...
        'Select files', ... 
        'MultiSelect', 'off');
 
if pathname == 0
    return
end

cd(pathname);
set(handles.FileNameBox, 'String', input_file);

LoadData(handles);

% --- MvdM utility functions
function UpdatePlot(handles)

global data;
global params;

axes(handles.axes2); cla; hold on;
if data.nEvents == 0, return; end

clear evt_data;
for iEvt = data.nEvents:-1:1
    
    temp_idx = data.all_tvec > data.evt_t(iEvt) + params.timeWindow(1) & data.all_tvec < data.evt_t(iEvt) + params.timeWindow(2);
    temp_data = data.all_lfp(temp_idx);
    
    if iEvt ~= data.nEvents
        % deal with fencepost errors
        if length(temp_data) > size(evt_data,2)
            temp_data = temp_data(1:size(evt_data,2));
        elseif length(temp_data) < size(evt_data,2)
            temp_data = cat(1,temp_data,nan(size(evt_data,2)-length(temp_data),1));
        end
    end
    
    evt_data(iEvt,:) = temp_data;
    
    plot(params.window_tvec,evt_data(iEvt,:)); hold on;
    
end
plot(params.window_tvec,nanmean(evt_data),'r','LineWidth',2);

set(gca,'XLim',params.timeWindow);
set(gca,'YLim',[str2double(get(handles.YLim1,'String')) str2double(get(handles.YLim2,'String'))]);

plot([0 0],ylim,'--k');


function LoadData(handles)

global data;
global params;

fn = get(handles.FileNameBox, 'String');

if ~isempty(fn)
    [data.all_lfp, data.all_tvec, ~] = load_open_ephys_data(fn);
    
    data.all_lfp = data.all_lfp(1:params.decimateFactor:end); % hack, but fast
    data.all_tvec = data.all_tvec(1:params.decimateFactor:end);
    
    [data.all_evt_id, data.all_evt_t, ~] = load_open_ephys_data('all_channels.events');
end

dt = median(diff(data.all_tvec));
params.window_tvec = params.timeWindow(1)+dt/2:dt:params.timeWindow(2);

% find event indices of sweep starts
data.sweep_evt_idx = find(data.all_evt_id == params.sweepChannel);
data.nSweeps = length(data.sweep_evt_idx);
set(handles.SweepBox,'String',data.nSweeps);

% find events for this sweep
FindEvents(handles);

% plot
UpdatePlot(handles);


function FindEvents(handles)

global data;
global params;

currentSweep = str2double(get(handles.SweepBox,'String'));
tstart = data.all_evt_t(data.sweep_evt_idx(currentSweep));

if currentSweep == data.nSweeps % last sweep
    tend = data.all_tvec(end);
else
    tend = data.all_evt_t(data.sweep_evt_idx(currentSweep+1));
end

data.evt_idx = find(data.all_evt_id == params.eventChannel);
data.evt_t = data.all_evt_t(data.evt_idx);
data.evt_t = data.evt_t(data.evt_t >= tstart & data.evt_t <= tend); 

data.nEvents = length(data.evt_t);
set(handles.EventBox,'String',data.nEvents);
    



function YLim2_Callback(hObject, eventdata, handles)
% hObject    handle to YLim2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of YLim2 as text
%        str2double(get(hObject,'String')) returns contents of YLim2 as a double
UpdatePlot(handles);

% --- Executes during object creation, after setting all properties.
function YLim2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to YLim2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function YLim1_Callback(hObject, eventdata, handles)
% hObject    handle to YLim1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of YLim1 as text
%        str2double(get(hObject,'String')) returns contents of YLim1 as a double
% find events for this sweep
UpdatePlot(handles);

% --- Executes during object creation, after setting all properties.
function YLim1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to YLim1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
