function varargout = proj(varargin)
% PROJ MATLAB code for proj.fig
%      PROJ, by itself, creates a new PROJ or raises the existing
%      singleton*.
%
%      H = PROJ returns the handle to a new PROJ or the handle to
%      the existing singleton*.
%
%      PROJ('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in PROJ.M with the given input arguments.
%
%      PROJ('Property','Value',...) creates a new PROJ or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before proj_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to proj_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help proj

% Last Modified by GUIDE v2.5 08-Apr-2013 15:09:08

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @proj_OpeningFcn, ...
                   'gui_OutputFcn',  @proj_OutputFcn, ...
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


% --- Executes just before proj is made visible.
function proj_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to proj (see VARARGIN)
clc;
World.points = generateworldpoints2();
setappdata(handles.figure1,'world',World);

Camera.f = 0.5;
Camera.ay = 480*Camera.f;
Camera.ax = 640*Camera.f;
Camera.u0 = 640/2;
Camera.v0 = 480/2;
Camera.Int = [Camera.ax 0 Camera.u0; 0 Camera.ay Camera.v0; 0 0 1];

Camera.camt = [0 0 0]';
Camera.thetax = 0;
Camera.thetay = 0;    
Camera.thetaz = 0;
Camera.Ext = [];
setappdata(handles.figure1,'camera',Camera);

InitKeyFrame1 = [];
InitKeyFrame2 = [];
setappdata(handles.figure1,'initkf1',InitKeyFrame1);
setappdata(handles.figure1,'initkf2',InitKeyFrame2);

Display(handles);






% Choose default command line output for proj
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes proj wait for user response (see UIRESUME)
% uiwait(handles.figure1);



% --- Outputs from this function are returned to the command line.
function varargout = proj_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in pushbutton_w.
function pushbutton_w_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_w (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
Camera = getappdata(handles.figure1,'camera');
ct = [0 0 1 0]';
ct = Camera.Ext\ct;
Camera.camt = Camera.camt + ct(1:3);
setappdata(handles.figure1,'camera',Camera);
Display(handles);


% --- Executes on button press in pushbutton_d.
function pushbutton_d_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_d (see GCBOfh = figure('Position',[250 250 350 35)
% eventdata  reserved - to be defined i4n a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
Camera = getappdata(handles.figure1,'camera');
ct = [1 0 0 0]';
ct = Camera.Ext\ct;
Camera.camt = Camera.camt + ct(1:3);
setappdata(handles.figure1,'camera',Camera);
Display(handles);

% --- Executes on button press in pushbutton_s.
function pushbutton_s_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_s (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
Camera = getappdata(handles.figure1,'camera');
ct = [0 0 -1 0]';
ct = Camera.Ext\ct;
Camera.camt = Camera.camt + ct(1:3);
setappdata(handles.figure1,'camera',Camera);
Display(handles);

% --- Executes on button press in pushbutton_a.
function pushbutton_a_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_a (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
Camera = getappdata(handles.figure1,'camera');
ct = [-1 0 0 0]';
ct = Camera.Ext\ct;
Camera.camt = Camera.camt + ct(1:3);
setappdata(handles.figure1,'camera',Camera);
Display(handles);
    
% --- Executes on button press in pushbutton_right.
function pushbutton_right_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_right (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handes and user data (see GUIDATA)
Camera = getappdata(handles.figure1,'camera');
Camera.thetay = Camera.thetay + 0.1;
guidata(handles.figure1,handles);
setappdata(handles.figure1,'camera',Camera);
Display(handles);


% --- Executes on button press in pucamerashbutton_left.
function pushbutton_left_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_left (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles anmydata = guidata(hObject);d user data (see GUIDATA)
Camera = getappdata(handles.figure1,'camera');
Camera.thetay = Camera.thetay - 0.1;
setappdata(handles.figure1,'camera',Camera);
Display(handles);


function Display(handles)
Camera = getappdata(handles.figure1,'camera');


Rx = [1 0 0; 0 cos(Camera.thetax) sin(Camera.thetax); 0 -sin(Camera.thetax) cos(Camera.thetax)];
Ry = [cos(Camera.thetay) 0 -sin(Camera.thetay); 0 1 0; sin(Camera.thetay) 0 cos(Camera.thetay)];
Rz = [cos(Camera.thetaz) sin(Camera.thetaz) 0; -sin(Camera.thetaz) cos(Camera.thetaz) 0; 0 0 1];
R = Rx*Ry*Rz;
t = -R*Camera.camt;
Camera.Ext = [R t; 0 0 0 1]; 
Pvanilla = [1 0 0 0; 0 1 0 0; 0 0 1 0];
Camera.P = Camera.Int*Pvanilla*Camera.Ext;
setappdata(handles.figure1,'camera',Camera);
World = getappdata(handles.figure1,'world');

CurrKeyFrame = MakeKeyFrame(Camera, World);
DisplayKeyFrame(CurrKeyFrame,handles.view3d);
setappdata(handles.figure1,'currkeyframe',CurrKeyFrame);


cla(handles.viewtopdown);
axes(handles.viewtopdown);
hold on;
plot(0, 0,'bx');
plot(0, 4,'wx');
plot(Camera.camt(1),Camera.camt(3),'gx');
Zaxis = (Camera.Ext)\[0 0 1 1]';
Xaxis = (Camera.Ext)\[1 0 0 1]';
plot(Zaxis(1),Zaxis(3),'bx');
plot(Xaxis(1),Xaxis(3),'rx');
hold off;


function [KeyFrame] = MakeKeyFrame(Camera, World)
kfimpointcount = 0;
KeyFrame.Camera = Camera;
KeyFrame.ImagePoints = struct('id',1,'location',[0 0 1]');
for i = 1:length(World.points)
    ImagePoint = ProjectPoint(Camera, World.points(i));
    if (~isempty(ImagePoint))
        kfimpointcount = kfimpointcount + 1;
        KeyFrame.ImagePoints(kfimpointcount) = ImagePoint;
    end
end


function [ImagePoint] = ProjectPoint(Camera, WorldPoint)
ImagePoint = [];
X = WorldPoint.location;
nX = Camera.Ext*X;
nX = nX ./ nX(4);
if (nX(3) > Camera.f)
    x = Camera.P*X;
    x = x./x(3);
    if (x(1) > 1 && x(1) < 640 && x(2) > 1 && x(2) < 480)
        ImagePoint.id = WorldPoint.id;
        ImagePoint.location = [x(1) x(2) 1]';
    end
end    


% --- Executes on button press in pushbutton_poke.
function pushbutton_poke_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_poke (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

InitKeyFrame1 = getappdata(handles.figure1,'initkf1');
InitKeyFrame2 = getappdata(handles.figure1,'initkf2');
CurrKeyFrame = getappdata(handles.figure1,'currkeyframe');

if (isempty(InitKeyFrame1))
    InitKeyFrame1 = CurrKeyFrame;
    DisplayKeyFrame(InitKeyFrame1, handles.viewkeyframe1);
    setappdata(handles.figure1,'initkf1',InitKeyFrame1);
else
    if (isempty(InitKeyFrame2))
        InitKeyFrame2 = CurrKeyFrame;
        DisplayKeyFrame(InitKeyFrame2, handles.viewkeyframe2);
        setappdata(handles.figure1,'initkf2',InitKeyFrame2);
    else
        Reproject(InitKeyFrame1,InitKeyFrame2);
    end

    
    
end

function Reproject(Keyframe1, Keyframe2)

kf1points = []; 
kf2points = [];

for i = 1:length(Keyframe1.ImagePoints)
    for j = 1:length(Keyframe2.ImagePoints)
        if (Keyframe1.ImagePoints(i).id == Keyframe2.ImagePoints(j).id)
            kf1points = [kf1points Keyframe1.ImagePoints(i).location];
            kf2points = [kf2points Keyframe2.ImagePoints(j).location];
        end
        
    end
end

X = [];
for i = 1:size(kf1points,2)
    X = [X linearreproject(kf1points(:,i),kf2points(:,i),Keyframe1.Camera.P,Keyframe2.Camera.P)];
end
display(X);


  

f = figure('Position',[100 100 640*2 480]);
set(gca,'Color',[0 0 0]);
axis([0 640*2 0 480]);
hold on;
plot([640 640],[0 480],'w-');
plot(kf1points(1,:), kf1points(2,:),'w+');
plot(640+kf2points(1,:), kf2points(2,:),'w+');
plot([kf1points(1,:); 640+kf2points(1,:)], [kf1points(2,:); kf2points(2,:)]);









function DisplayKeyFrame(KeyFrame, AxesHandle)
cla(AxesHandle);
axes(AxesHandle);
hold on;
for i = 1:length(KeyFrame.ImagePoints)
    plot(KeyFrame.ImagePoints(i).location(1), KeyFrame.ImagePoints(i).location(2),'w');
end

%set(AxesHandle,'Color',[0 0 0]);
%set(AxesHandle,'XLim',[0 640]);
%set(AxesHandle,'YLim',[0 480]);


% --- Executes on button press in pushbutton_path.
function pushbutton_path_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_path (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
Camera = getappdata(handles.figure1,'camera');
for theta = 0:0.1:2*pi
    clc
    ct = [9*cos(theta) 0 9*sin(theta) 1]';
    display(ct);
    Camera.camt = ct(1:3);
    Camera.thetay = -theta;
    setappdata(handles.figure1,'camera',Camera);
    Display(handles);
end

