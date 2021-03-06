function outputParams = simpleBlendGUI(source, target)
%------------------------------------------------
% Set up the GUI
% Written by: Aria Pezeshk, March 2015
%------------------------------------------------
h_tmp = figure('Visible','Off'); %this is a temp figure, just so that when 'Blend' callback closes it, gui can be closed too and outputs returned
h_gui = figure('units','normalized','outerposition',[0 0 1 0.5]);
set(h_gui, 'Units', 'Pixels'); %revert to default pixels; somehow imcrop, imrect, impoint, ... don't work unless units is pixels!
handles.h_gui = h_gui;
handles.axesSource = axes('Parent', h_gui, ... % Axes for displaying source image
        'Units', 'normalized', ...
        'Tag', 'axesSource', ...        
        'Position',[0.05 0.05 0.4 0.8]);         
handles.axesTarget = axes('Parent', h_gui, ... % Axes for displaying source image
        'Units', 'normalized', ...                 
        'Tag', 'axesTarget', ...
        'Position',[0.55 0.05 0.4 0.8]);         
handles.pushButtonCropROI = uicontrol('Parent', h_gui, ...
        'Style', 'pushbutton', ...
        'String', 'Crop ROI', ...
        'Units', 'normalized', ...
        'Position', [0.47 0.5 0.05 0.05],...
        'Tag', 'pushButtonCropROI', ...
        'Callback', @CropROI_Callback);  
handles.pushButtonRoughSegment = uicontrol('Parent', h_gui, ...
        'Style', 'pushbutton', ...
        'String', 'Rough Segment', ...
        'Units', 'normalized', ...
        'Position', [0.47 0.4 0.05 0.05],...
        'Tag', 'pushButtonRoughSegment', ...
        'Callback', @RoughSegment_Callback);
handles.pushButtonSelectTarget = uicontrol('Parent', h_gui, ...
        'Style', 'pushbutton', ...
        'String', 'Select Target', ...
        'Units', 'normalized', ...
        'Position', [0.47 0.3 0.05 0.05],...
        'Tag', 'pushButtonSelectTarget', ...
        'Callback', @SelectTarget_Callback);
handles.pushButtonBlend = uicontrol('Parent', h_gui, ...
        'Style', 'pushbutton', ...
        'String', 'Blend', ...
        'Units', 'normalized', ...
        'Position', [0.47 0.2 0.05 0.05],...
        'Tag', 'pushButtonBlend', ...
        'Callback', @Blend_Callback);


handles.sourceImage = imshow(source, 'Parent', handles.axesSource, 'DisplayRange',[]);  
handles.srcNumRows = size(source,1); handles.srcNumCols = size(source, 2);
handles.targetImage = imshow(target, 'Parent', handles.axesTarget, 'DisplayRange',[]);  
handles.rectPos = []; %initialize
handles.pointPos = []; %initialize
handles.segMaskFullSize = []; %initialize
handles.hPoint = []; %initialize
handles.h_tmp = h_tmp;
guidata(h_gui, handles);

waitfor(h_tmp);
handles = guidata(h_gui);
outputParams.rectPos = handles.rectPos;
outputParams.pointPos = handles.pointPos;
outputParams.segMaskFullSize = handles.segMaskFullSize;

delete(h_gui);

%------------------------------------------------
function CropROI_Callback(hObject, callbackdata)
handles = guidata(hObject);
zoom off; % in case zoom mode is on, turn it off so that interactive mode is exited and following can be executed without error
datacursormode off;% in case datacursormode is on, turn it off so that interactive mode is exited and following can be executed without error
pan off;

% hRect = imrect(handles.axesSource,[]);
% API = iptgetapi(hRect);
% handles.rectPos = round(API.getPosition());
% delete(hRect);
[tmp, rectPos] = imcrop(handles.sourceImage);
handles.rectPos = round(rectPos);
guidata(hObject, handles);

%------------------------------------------------
function RoughSegment_Callback(hObject, callbackdata)
handles = guidata(hObject);
hPolyLesion = impoly(handles.axesSource,[]);
API = iptgetapi(hPolyLesion);
polyVertices= API.getPosition();
delete(hPolyLesion);
handles.segMaskFullSize = poly2mask(polyVertices(:,1), polyVertices(:,2), handles.srcNumRows, handles.srcNumCols); %this gives a mask image with same size as source
guidata(hObject, handles);


%------------------------------------------------
function SelectTarget_Callback(hObject, callbackdata)
handles = guidata(hObject);
zoom off; % in case zoom mode is on, turn it off so that interactive mode is exited and following can be executed without error
datacursormode off;% in case datacursormode is on, turn it off so that interactive mode is exited and following can be executed without error
pan off;

hPointOld = handles.hPoint; % check if there is still a previous point on the figure, and if so delete it
if ishandle(hPointOld)
        delete(hPointOld);
end
hPoint = impoint(handles.axesTarget,[]);
API = iptgetapi(hPoint);
API.setColor('r');
handles.hPoint = hPoint;
API = iptgetapi(hPoint);
handles.pointPos= round(API.getPosition());
guidata(hObject,handles);


%------------------------------------------------
function Blend_Callback(hObject, callbackdata)
handles = guidata(hObject);
delete(handles.h_tmp);
