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

% Last Modified by GUIDE v2.5 18-Apr-2013 14:38:03

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
ay = 480*Camera.f;
ax = 640*Camera.f;
u0 = 640/2;
v0 = 480/2;
Camera.K = [ax 0 u0; 0 ay v0; 0 0 1];

Camera.camt = [0 0 0]';
Camera.thetax = 0;
Camera.thetay = 0;    
Camera.thetaz = 0;
Camera.R = eye(3,3);
Camera.t = zeros(3,1);
Camera.E = [Camera.R Camera.t; 0 0 0 1];
Camera.P = Camera.K*Camera.E(1:3,:);
setappdata(handles.figure1,'camera',Camera);
EstCamera = Camera;
setappdata(handles.figure1,'estcamera',EstCamera);

InitKeyFrame1 = [];
InitKeyFrame2 = [];
setappdata(handles.figure1,'initkf1',InitKeyFrame1);
setappdata(handles.figure1,'initkf2',InitKeyFrame2);


EstWorld.points(1).location = [0 0 0 1]';
EstWorld.points(1).id = -1;
% EstWorld.points = World.points;
setappdata(handles.figure1,'estworld',EstWorld);

UpdateTick(handles)







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
ct = Camera.E\ct;
Camera.camt = Camera.camt + ct(1:3);
setappdata(handles.figure1,'camera',Camera);
UpdateTick(handles);


% --- Executes on button press in pushbutton_d.
function pushbutton_d_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_d (see GCBOfh = figure('Position',[250 250 350 35)
% eventdata  reserved - to be defined i4n a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
Camera = getappdata(handles.figure1,'camera');
ct = [1 0 0 0]';
ct = Camera.E\ct;
Camera.camt = Camera.camt + ct(1:3);
setappdata(handles.figure1,'camera',Camera);
UpdateTick(handles);

% --- Executes on button press in pushbutton_s.
function pushbutton_s_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_s (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
Camera = getappdata(handles.figure1,'camera');
ct = [0 0 -1 0]';
ct = Camera.E\ct;
Camera.camt = Camera.camt + ct(1:3);
setappdata(handles.figure1,'camera',Camera);
UpdateTick(handles);

% --- Executes on button press in pushbutton_a.
function pushbutton_a_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_a (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
Camera = getappdata(handles.figure1,'camera');
ct = [-1 0 0 0]';
ct = Camera.E\ct;
Camera.camt = Camera.camt + ct(1:3);
setappdata(handles.figure1,'camera',Camera);
UpdateTick(handles);

% --- Executes on button press in pushbutton_right.
function pushbutton_right_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_right (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handes and user data (see GUIDATA)
Camera = getappdata(handles.figure1,'camera');
Camera.thetay = Camera.thetay + 0.1;
guidata(handles.figure1,handles);
setappdata(handles.figure1,'camera',Camera);
UpdateTick(handles);



% --- Executes on button press in pucamerashbutton_left.
function pushbutton_left_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_left (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles anmydata = guidata(hObject);d user data (see GUIDATA)
Camera = getappdata(handles.figure1,'camera');
Camera.thetay = Camera.thetay - 0.1;
setappdata(handles.figure1,'camera',Camera);
UpdateTick(handles);

function outputCamera = RfromEuler(Camera)
outputCamera = Camera;
Rx = [1 0 0; 0 cos(Camera.thetax) sin(Camera.thetax); 0 -sin(Camera.thetax) cos(Camera.thetax)];
Ry = [cos(Camera.thetay) 0 -sin(Camera.thetay); 0 1 0; sin(Camera.thetay) 0 cos(Camera.thetay)];
Rz = [cos(Camera.thetaz) sin(Camera.thetaz) 0; -sin(Camera.thetaz) cos(Camera.thetaz) 0; 0 0 1];
outputCamera.R = Rx*Ry*Rz;
outputCamera.t = -outputCamera.R*Camera.camt;
outputCamera.E = [outputCamera.R outputCamera.t; 0 0 0 1];

function DisplayTopDown(Camera, viewhandle)
cla(viewhandle);
axes(viewhandle);
hold on;
plot(0, 0,'bx');
Yaxis = (Camera.E)\[0 1 0 1]';
Zaxis = (Camera.E)\[0 0 1 1]';
Xaxis = (Camera.E)\[1 0 0 1]';
plot(Zaxis(1),Zaxis(3),'bx');
plot(Yaxis(1),Yaxis(3),'gx');
plot(Xaxis(1),Xaxis(3),'rx');
hold off;



function [ImagePoints] = MakeImage(Camera, World)
kfimpointcount = 0;
ImagePoints = struct('id',1,'location',[0 0 1]','X',[0 0 0 1]');
for i = 1:length(World.points)
    ImagePoint = ProjectPoint(Camera, World.points(i));
    if (~isempty(ImagePoint))
        kfimpointcount = kfimpointcount + 1;
        ImagePoints(kfimpointcount) = ImagePoint;
    end
end


function [ImagePoint] = ProjectPoint(Camera, WorldPoint)
ImagePoint = [];
X = WorldPoint.location;
nX = Camera.E*X;
nX = nX ./ nX(4);
if (nX(3) > Camera.f)
    x = Camera.K*Camera.E(1:3,:)*X;
    x = x./x(3);
    if (x(1) > 1 && x(1) < 640 && x(2) > 1 && x(2) < 480)
        ImagePoint.id = WorldPoint.id;
        ImagePoint.location = [x(1) x(2) 1]';
        ImagePoint.X = X;
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
World = getappdata(handles.figure1,'world');

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
        points = Reproject(InitKeyFrame1,InitKeyFrame2);
        EstWorld.points = points;
        error = calculateworlderror(World,EstWorld);
        display(error);
    end

    
    
end

function Ext = CalculateExt(Keyframe1, Keyframe2,K, F1)

kf1points = []; 
kf2points = [];
ids = [];
ids2 = [];

for i = 1:length(Keyframe1.ImagePoints)
    for j = 1:length(Keyframe2.ImagePoints)
        if (Keyframe1.ImagePoints(i).id == Keyframe2.ImagePoints(j).id)
            kf1points = [kf1points Keyframe1.ImagePoints(i).location];
            kf2points = [kf2points Keyframe2.ImagePoints(j).location];
            ids = [ids Keyframe1.ImagePoints(i).id];
            ids2 = [ids2 Keyframe2.ImagePoints(j).id];
            
        end
        
    end
end


F = fundmatrix(kf1points,kf2points);
% F = vgg_F_from_7pts_2img(kf1points(:,1:7),kf2points(:,1:7));
display(F);
E = K'*F*K;
t = null(E');
display(t);
t1 = t;
t2 = -t;




[U, S, V] = svd(E);

W = [0 -1 0; 1 0 0; 0 0 1];

R1 = U*W*V';
R2 = U*W'*V';

% display(R1);
% display(R2);

P = K*[eye(3,3) zeros(3,1)];
P1 = K*[R1 t1];
P2 = K*[R1 t2];
P3 = K*[R2 t1];
P4 = K*[R2 t2];




error1 = 0;
error2 = 0;
error3 = 0;
error4 = 0;

for i = 1:size(kf1points,2)
    X1 = linearreproject(kf1points(:,i),kf2points(:,i),P,P1);
    X2 = linearreproject(kf1points(:,i),kf2points(:,i),P,P2);
    X3 = linearreproject(kf1points(:,i),kf2points(:,i),P,P3);
    X4 = linearreproject(kf1points(:,i),kf2points(:,i),P,P4);
    
    if (X1(3)>0)
        X1 = [R1 t1; 0 0 0 1]*X1;
        error1 = error1 + (X1(3)<0);
    else
        error1 = error1 + 1;
    end
    
    if (X2(3)>0)
        X2 = [R1 t2; 0 0 0 1]*X2;
        error2 = error2 + (X2(3)<0);
    else
        error2 = error2 + 1;
    end
    
    if (X3(3)>0)
        X3 = [R2 t1; 0 0 0 1]*X3;
        error3 = error3 + (X3(3)<0);
    else
        error3 = error3 + 1;
    end
    
    if (X4(3)>0)
        X4 = [R2 t2; 0 0 0 1]*X4;
        error4 = error4 + (X4(3)<0);
    else
        error4 = error4 + 1;
    end    
   
    
end

display(error1);
display(error2);
display(error3);
display(error4);

if error1 == 0
    R = R1;
    t = t1;
end

if error2 == 0
    R = R1;
    t = t2;
end

if error3 == 0
    R = R2;
    t = t1;
end

if error4 == 0
    R = R2;
    t = t2;
end

scale = abs(2/t(1));
t = t*scale;



display(R);
display(t);

Ext = [R t; 0 0 0 1];





















function outpoints = Reproject(Keyframe1, Keyframe2, P1, P2)

kf1points = []; 
kf2points = [];
ids = [];

for i = 1:length(Keyframe1.ImagePoints)
    for j = 1:length(Keyframe2.ImagePoints)
        if (Keyframe1.ImagePoints(i).id == Keyframe2.ImagePoints(j).id)
            kf1points = [kf1points Keyframe1.ImagePoints(i).location];
            kf2points = [kf2points Keyframe2.ImagePoints(j).location];
            ids = [ids Keyframe1.ImagePoints(i).id];
            
        end
        
    end
end

X = [];

for i = 1:size(kf1points,2)
    newpoint = linearreproject(kf1points(:,i),kf2points(:,i),P1,P2);
    X = [X newpoint];
    outpoints(i).location = newpoint;
    outpoints(i).id = ids(i);
end

function DisplayImage(ImagePoints, AxesHandle)
cla(AxesHandle);
axes(AxesHandle);
hold on;
for i = 1:length(ImagePoints)
    plot(ImagePoints(i).location(1), ImagePoints(i).location(2),'w');
end


% --- Executes on button press in pushbutton_path.
function pushbutton_path_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_path (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

Camera = getappdata(handles.figure1,'camera');
EstCamera = getappdata(handles.figure1,'estcamera');
errcount = 0;
err = [];
err2 = [];
description = 'Relaxed error';
for theta = 0:0.1:2*pi
    EstCamera = getappdata(handles.figure1,'estcamera');
    errcount = errcount + 1;
    err2(errcount) = norm(Camera.E - EstCamera.E);
    clc
    ct = [18*cos(theta)-16 0 18*sin(theta)]';
    Camera.camt = ct;
    Camera.thetay = -theta;
    Camera = RfromEuler(Camera);

    setappdata(handles.figure1,'camera',Camera);
    UpdateTick(handles);
    
    err(errcount) = EstimateCamera(handles);
    save results9 err err2 description theta EstCamera Camera
end


% --- Executes on button press in pushbutton_init1.
function pushbutton_init1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_init1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
Camera = getappdata(handles.figure1,'camera');
ct = [0 0 0]';
Camera.camt = ct;
Camera.thetay = 0;
setappdata(handles.figure1,'camera',Camera);
UpdateTick(handles);

CurrKeyFrame = getappdata(handles.figure1,'currkeyframe');
InitKeyFrame1 = CurrKeyFrame;
InitKeyFrame1.Camera = Camera;
DisplayImage(InitKeyFrame1.ImagePoints, handles.viewkeyframe1);
setappdata(handles.figure1,'initkf1',InitKeyFrame1);
setappdata(handles.figure1,'camera',Camera);



% --- Executes on button press in pushbutton_init2.
function pushbutton_init2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_init2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
Camera = getappdata(handles.figure1,'camera');
ct = [2 0 0]';
Camera.camt = ct;
Camera.thetay = -0.1;
setappdata(handles.figure1,'camera',Camera);
UpdateTick(handles);

CurrKeyFrame = getappdata(handles.figure1,'currkeyframe');
InitKeyFrame2 = CurrKeyFrame;
DisplayImage(InitKeyFrame2.ImagePoints, handles.viewkeyframe2);
setappdata(handles.figure1,'initkf2',InitKeyFrame2);



% --- Executes on button press in pushbutton_reproject.
function pushbutton_reproject_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_reproject (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
InitKeyFrame1 = getappdata(handles.figure1,'initkf1');
InitKeyFrame2 = getappdata(handles.figure1,'initkf2');
World = getappdata(handles.figure1,'world');
Camera = getappdata(handles.figure1,'camera');
EstCamera = getappdata(handles.figure1,'estcamera');

display(Camera.P);
tx = [0 0 2; 0 0 0; -2 0 0];
E1 = 0;
display(Camera.R);
EstCamera.E = CalculateExt(InitKeyFrame1, InitKeyFrame2, Camera.K,E1);
EstCamera.P = EstCamera.K*EstCamera.E(1:3,:);
P1 = EstCamera.K*[eye(3,3) [0;0;0]];

InitKeyFrame2.Camera = EstCamera;


newpoints = Reproject(InitKeyFrame1, InitKeyFrame2,P1,EstCamera.P);
EstWorld.points = newpoints;
setappdata(handles.figure1,'estworld',EstWorld);
setappdata(handles.figure1,'estcamera',EstCamera);
UpdateTick(handles);
error = calculateworlderror(World,EstWorld);
set(handles.text_worlderror,'String',['World Error: ' num2str(error)]);

function UpdateTick(handles)
World = getappdata(handles.figure1,'world');
EstWorld = getappdata(handles.figure1,'estworld');
Camera = getappdata(handles.figure1,'camera');
EstCamera = getappdata(handles.figure1,'estcamera');
Camera = RfromEuler(Camera);


%Display the camera external matrices and the error 
E1 = Camera.E;
E1 = round(E1*1000)/1000;
set(handles.text_camext,'String',['Camera E: ' mat2str(E1)]);

E2 = EstCamera.E;
E2 = round(E2*1000)/1000;
set(handles.text_estcamext,'String',['EstCamera E: ' mat2str(E2)]);

Cam_Error = norm(E2 - E1);
set(handles.text_cameraerror,'String',['Camera Error: ' num2str(Cam_Error)]);


setappdata(handles.figure1,'camera',Camera);
setappdata(handles.figure1,'estcamera',EstCamera);


CurrKeyFrame.ImagePoints = MakeImage(Camera, World);
DisplayImage(CurrKeyFrame.ImagePoints, handles.view3d);



setappdata(handles.figure1,'currkeyframe',CurrKeyFrame);
DisplayTopDown(Camera,handles.viewtopdown);

EstImagePoints = MakeImage(EstCamera, EstWorld);
DisplayImage(EstImagePoints, handles.view3dest);
% setappdata(handles.figure1,'estcurrkeyframe',EstCurrKeyFrame);
DisplayTopDown(EstCamera,handles.viewtopdownest);



% --- Executes on button press in pushbutton_estcam.
function pushbutton_estcam_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_estcam (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
EstimateCamera(handles);





% --- Executes on button press in pushbutton_update.
function pushbutton_update_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_update (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
UpdateTick(handles);

function error = EstimateCamera(handles)
World = getappdata(handles.figure1,'world');
EstWorld = getappdata(handles.figure1,'estworld');
Camera = getappdata(handles.figure1,'camera');
EstCamera = getappdata(handles.figure1,'estcamera');

CurrKeyFrame.ImagePoints = MakeImage(Camera, World);



Ein = EstCamera.E;

K = EstCamera.K;
XX = [];
detection = [];
npoints = 50;
error = 10;


[X x] = FindMatches(EstWorld, CurrKeyFrame,npoints);
if (size(X,2) > 10)
    display('Estimating pose...');
    [muout error] = estimatepose(Ein,K,X,x,300);
    EstCamera.E = expmap(muout)*Ein;
    EstCamera.mu = muout;
    EstCamera.P = K*EstCamera.E(1:3,:);
    display('Done estimating!');
else
    display('Not enough points to estimate pose');
end

setappdata(handles.figure1,'estcamera',EstCamera);
UpdateTick(handles);


function [X x] = FindMatches(World, KeyFrame,npoints)
X = [];
x = [];

for i = 1:size(KeyFrame.ImagePoints,2)
    input(1:3,i) = KeyFrame.ImagePoints(i).location;
    input(4,i) = KeyFrame.ImagePoints(i).id;
end

input = input(:,randperm(size(input,2)));

for i = 1:npoints
    for j = 1:size(World.points,2)
        if (World.points(j).id == input(4,i))
            X = [X World.points(j).location];
            x = [x input(1:3,i)];
        end
    end
end


% --- Executes on button press in pushbutton_addkeyframe.
function pushbutton_addkeyframe_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_addkeyframe (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


