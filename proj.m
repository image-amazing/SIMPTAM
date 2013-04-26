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

% Last Modified by GUIDE v2.5 24-Apr-2013 13:45:40

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
Map.points = generateworldpoints2();

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

World.Camera = Camera;
World.Map = Map;
PTAM.Camera = Camera;


PTAM.Map = [];
PTAM.Map.points(1).location = [0 0 0 0]';
PTAM.Map.points(1).id = -1;
PTAM.kfcount = 0;

% PTAM.KeyFrames = struct([]);

setappdata(handles.figure1,'world',World);
setappdata(handles.figure1,'ptam',PTAM);

theta = 0;
setappdata(handles.figure1,'theta',theta);

InitKeyFrame1 = [];
InitKeyFrame2 = [];


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
World = getappdata(handles.figure1,'world');
ct = [0 0 1 0]';
ct = World.Camera.E\ct;
World.Camera.camt = World.Camera.camt + ct(1:3);
setappdata(handles.figure1,'world',World);
UpdateTick(handles);


% --- Executes on button press in pushbutton_d.
function pushbutton_d_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_d (see GCBOfh = figure('Position',[250 250 350 35)
% eventdata  reserved - to be defined i4n a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
World = getappdata(handles.figure1,'world');
ct = [1 0 0 0]';
ct =  World.Camera.E\ct;
World.Camera.camt =  World.Camera.camt + ct(1:3);
setappdata(handles.figure1,'world',World);
UpdateTick(handles);

% --- Executes on button press in pushbutton_s.
function pushbutton_s_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_s (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
World = getappdata(handles.figure1,'world');
ct = [0 0 -1 0]';
ct = World.Camera.E\ct;
World.Camera.camt = World.Camera.camt + ct(1:3);
setappdata(handles.figure1,'world',World);
UpdateTick(handles);

% --- Executes on button press in pushbutton_a.
function pushbutton_a_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_a (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
World = getappdata(handles.figure1,'world');
ct = [-1 0 0 0]';
ct = World.Camera.E\ct;
World.Camera.camt = World.Camera.camt + ct(1:3);
setappdata(handles.figure1,'world',World);
UpdateTick(handles);

% --- Executes on button press in pushbutton_right.
function pushbutton_right_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_right (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handes and user data (see GUIDATA)
World = getappdata(handles.figure1,'world');
World.Camera.thetay = World.Camera.thetay + 0.1;
guidata(handles.figure1,handles);
setappdata(handles.figure1,'world',World);
UpdateTick(handles);



% --- Executes on button press in pucamerashbutton_left.
function pushbutton_left_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_left (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles anmydata = guidata(hObject);d user data (see GUIDATA)
World = getappdata(handles.figure1,'world');
World.Camera.thetay = World.Camera.thetay - 0.1;
setappdata(handles.figure1,'world',World);
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

function DisplayKeyFrames(KeyFrames, viewhandle)
axes(viewhandle);
hold on;

for i = 1:size(KeyFrames,2)
    Yaxis = (KeyFrames(i).Camera.E)\[0 1 0 1]';
    Zaxis = (KeyFrames(i).Camera.E)\[0 0 1 1]';
    Xaxis = (KeyFrames(i).Camera.E)\[1 0 0 1]';
    plot(Xaxis(1),Xaxis(3),'x','Color',[0.5 0 0]);
    plot(Yaxis(1),Yaxis(3),'x','Color',[0 0.5 0]);
    plot(Zaxis(1),Zaxis(3),'x','Color',[0 0 0.5]);
end
hold off;



function [ImagePoints] = MakeImage(Camera, World)
kfimpointcount = 0;
ImagePoints = struct('id',1,'location',[0 0 1]','location2',[0 0 1]' ,'X',[0 0 0 1]');
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
    
    x2 = double(Camera.E(1:3,:)*X);
    x2 = x2./x2(3);
    
    if (x(1) > 1 && x(1) < 640 && x(2) > 1 && x(2) < 480)
        ImagePoint.id = WorldPoint.id;
        ImagePoint.location = [x(1) x(2) 1]';
        ImagePoint.location2 = [double(x2(1)) double(x2(2)) 1]';
        ImagePoint.X = X;
    end
end


% --- Executes on button press in pushbutton_poke.
function pushbutton_poke_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_poke (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


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

function Ext = CalculateExt(Keyframe1, Keyframe2,K)

kf1points = [];
kf2points = [];
ids = [];
ids2 = [];

for i = 1:length(Keyframe1.ImagePoints)
    for j = 1:length(Keyframe2.ImagePoints)
        if (Keyframe1.ImagePoints(i).id == Keyframe2.ImagePoints(j).id)
            kf1p = Keyframe1.ImagePoints(i).location;
            kf1p = kf1p / kf1p(3);
            kf2p = Keyframe2.ImagePoints(j).location;
            kf2p = kf2p / kf2p(3);
            kf1points = [kf1points kf1p];
            kf2points = [kf2points kf2p];
            ids = [ids Keyframe1.ImagePoints(i).id];
            ids2 = [ids2 Keyframe2.ImagePoints(j).id];
            
        end
        
    end
end

R = [];
t = [];

while isempty(R)
    
    
    perm = randperm(size(kf1points,2));
    kf1points = kf1points(:,perm);
    kf2points = kf2points(:,perm);
    % kf1points = kf1points(:,1:12);
    % kf2points = kf2points(:,1:12);
    
    F = fundmatrix(kf1points,kf2points);
    % F = vgg_F_from_7pts_2img(kf1points(:,1:7),kf2points(:,1:7));
    display(F);
    E = K'*F*K;
    t = null(E');
    display(t);
    t1 = t;
    t2 = -t;
    
    ferror = 0;
    for i = 1:size(kf1points,2)
        ferror = ferror + abs(kf2points(:,i)'*F*kf1points(:,i));
    end
    display(ferror);
    
    
    
    
    [U, S, V] = svd(E);
    
    W = [0 -1 0; 1 0 0; 0 0 1];
    
    R1 = U*W*V';
    R2 = U*W'*V';
    
    % display(R1);
    % display(R2);
    
    P = [eye(3,3) zeros(3,1)];
    P1 = [R1 t1];
    P2 = [R1 t2];
    P3 = [R2 t1];
    P4 = [R2 t2];
    
    
    
    
    error1 = 0;
    error2 = 0;
    error3 = 0;
    error4 = 0;
    
    for i = 1:size(kf1points,2)
        X1 = linearreproject(K\kf1points(:,i),K\kf2points(:,i),P,P1);
        X2 = linearreproject(K\kf1points(:,i),K\kf2points(:,i),P,P2);
        X3 = linearreproject(K\kf1points(:,i),K\kf2points(:,i),P,P3);
        X4 = linearreproject(K\kf1points(:,i),K\kf2points(:,i),P,P4);
        
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
    
end

scale = abs(2/t(1));
t = t*scale;



display(R);
display(t);

R = single(R);
t = single(t);


Ext = [R t; 0 0 0 1];


function outpoints = Reproject(Keyframe1, Keyframe2, E1, E2,K)

kf1points = [];
kf2points = [];
ids = [];
X = [];

for i = 1:length(Keyframe1.ImagePoints)
    for j = 1:length(Keyframe2.ImagePoints)
        if (Keyframe1.ImagePoints(i).id == Keyframe2.ImagePoints(j).id)
            
            kf1p = K\Keyframe1.ImagePoints(i).location;
            kf1p = kf1p / kf1p(3);
            kf1p(1:2) = kf1p(1:2) + 0.0005*rand(2,1);
            kf2p = K\Keyframe2.ImagePoints(j).location;
            kf2p = kf2p / kf2p(3);
            kf1p(1:2) = kf1p(1:2) + 0.0005*rand(2,1);
            
            
            kf1points = [kf1points kf1p];
            kf2points = [kf2points kf2p];
            ids = [ids Keyframe1.ImagePoints(i).id];
            %             X = [X inpoints(Keyframe1.ImagePoints(i).id).location];
            
        end
        
    end
end

outpoints = [];
XX = [];
for i = 1:size(kf1points,2)
    newpoint = linearreproject(kf1points(:,i),kf2points(:,i),E1,E2);
    if ~isempty(newpoint)
        XX = [XX newpoint];
        outpoints(i).location = newpoint;
        outpoints(i).id = ids(i);
    end
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


keycount = 0;

for theta = 0.1:0.1:pi/2
    World = getappdata(handles.figure1,'world');
    keycount = keycount + 1;
    UpdateTick(handles);
    ct = [18*cos(theta)-16 0 18*sin(theta)]';
    World.Camera.camt = ct;
    World.Camera.thetay = -theta;
    World.Camera = RfromEuler( World.Camera);
    setappdata(handles.figure1,'world',World);
    setappdata(handles.figure1,'theta',theta);
    
    EstimateCamera(handles);
    UpdateTick(handles);
    if keycount > 0
        AddKeyFrame(handles);
        keycount = 0;
    end
    
end


% --- Executes on button press in pushbutton_init1.
function pushbutton_init1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_init1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
World = getappdata(handles.figure1,'world');
PTAM = getappdata(handles.figure1,'ptam');
ct = [0 0 0]';
World.Camera.camt = ct;
World.Camera.thetay = 0;
setappdata(handles.figure1,'world',World);
UpdateTick(handles);
World = getappdata(handles.figure1,'world');

CurrKeyFrame = getappdata(handles.figure1,'currkeyframe');

PTAM.KeyFrames(1) = CurrKeyFrame;
CurrKeyFrame.Camera = World.Camera;
World.KeyFrames(1) = CurrKeyFrame;
PTAM.kfcount = 1;
DisplayImage(PTAM.KeyFrames(1).ImagePoints, handles.viewkeyframe1);
setappdata(handles.figure1,'ptam',PTAM);
setappdata(handles.figure1,'world',World);

% --- Executes on button press in pushbutton_init2.
function pushbutton_init2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_init2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
World = getappdata(handles.figure1,'world');
PTAM = getappdata(handles.figure1,'ptam');
ct = [2 0 0]';
World.Camera.camt = ct;
World.Camera.thetay = -0.1;
setappdata(handles.figure1,'world',World);
UpdateTick(handles);
World = getappdata(handles.figure1,'world');

CurrKeyFrame = getappdata(handles.figure1,'currkeyframe');

PTAM.KeyFrames(2) = CurrKeyFrame;
CurrKeyFrame.Camera = World.Camera;
World.KeyFrames(2) = CurrKeyFrame;
PTAM.kfcount = 2;
DisplayImage(PTAM.KeyFrames(2).ImagePoints, handles.viewkeyframe2);
setappdata(handles.figure1,'ptam',PTAM);
setappdata(handles.figure1,'world',World);

% --- Executes on button press in pushbutton_reproject.
function pushbutton_reproject_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_reproject (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
World = getappdata(handles.figure1,'world');
PTAM = getappdata(handles.figure1,'ptam');

K = PTAM.Camera.K;
display(World.Camera.R);
PTAM.Camera = World.Camera;
PTAM.Camera.E = CalculateExt(PTAM.KeyFrames(1), PTAM.KeyFrames(2), World.Camera.K);

E1 = [eye(3,3) [0;0;0]];
E2 = World.Camera.E(1:3,:);



PTAM.KeyFrames(2).Camera = PTAM.Camera;   


PTAM.Map.points = Reproject(PTAM.KeyFrames(1), PTAM.KeyFrames(2),E1,E2,K);

setappdata(handles.figure1,'ptam',PTAM);
UpdateTick(handles);


function UpdateTick(handles)
World = getappdata(handles.figure1,'world');
PTAM = getappdata(handles.figure1,'ptam');
World.Camera = RfromEuler(World.Camera);


%Display the camera external matrices and the error
E1 = World.Camera.E;
E1 = round(E1*1000)/1000;
set(handles.text_camext,'String',['World Camera E: ' mat2str(E1)]);

E2 = PTAM.Camera.E;
E2 = round(E2*1000)/1000;
set(handles.text_estcamext,'String',['PTAM Camera E: ' mat2str(E2)]);

camerror = norm(E2 - E1);
set(handles.text_cameraerror,'String',['Camera Error: ' num2str(camerror)]);


if PTAM.kfcount > 1
    totalcamerror = 0;
    for i = 1:PTAM.kfcount
        E1 = World.KeyFrames(i).Camera.E;
        E2 = PTAM.KeyFrames(i).Camera.E;
        totalcamerror = totalcamerror + norm(E1-E2);
    end
    set(handles.text_totalcamerror,'String',['Total Camera Error: ' num2str(totalcamerror)]);
    set(handles.text_numkeyframes,'String',['Number of KeyFrames: ' num2str(PTAM.kfcount )]);
    set(handles.text_averagecamerror,'String',['Average Camera Error: ' num2str(totalcamerror/PTAM.kfcount )]);
end



[error count] = calculateworlderror(World.Map,PTAM.Map);
set(handles.text_totalmaperror,'String',['Total Map Error: ' num2str(error)]);
set(handles.text_mappoints,'String',['Number of Map Points: ' num2str(count)]);
set(handles.text_averagemaperror,'String',['Average Map Error: ' num2str(error/count)]);


setappdata(handles.figure1,'world',World);
setappdata(handles.figure1,'ptam',PTAM);


CurrKeyFrame.ImagePoints = MakeImage(World.Camera, World.Map);
CurrKeyFrame.Camera = PTAM.Camera;
DisplayImage(CurrKeyFrame.ImagePoints, handles.view3d);



setappdata(handles.figure1,'currkeyframe',CurrKeyFrame);
DisplayTopDown(World.Camera,handles.viewtopdown);

EstImagePoints = MakeImage(PTAM.Camera, PTAM.Map);
DisplayImage(EstImagePoints, handles.view3dest);
% setappdata(handles.figure1,'estcurrkeyframe',EstCurrKeyFrame);
DisplayTopDown(PTAM.Camera,handles.viewtopdownest);

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
PTAM = getappdata(handles.figure1,'ptam');

CurrKeyFrame.ImagePoints = MakeImage(World.Camera, World.Map);



Ein = PTAM.Camera.E;

K = PTAM.Camera.K;
XX = [];
detection = [];
npoints = 70;
error = 10;


[X x] = FindMatches(PTAM.Map, CurrKeyFrame,npoints);
if (size(X,2) > 10)
    display('Estimating pose...');
    [muout error] = estimatepose(Ein,K,X,x,300);
    PTAM.Camera.E = expmap(muout)*Ein;
%     PTAM.Camera.E = World.Camera.E;
    PTAM.Camera.mu = muout;
    PTAM.Camera.P = K*PTAM.Camera.E(1:3,:);
    display('Done estimating!');
else
    display('Not enough points to estimate pose');
end

setappdata(handles.figure1,'ptam',PTAM);

function [X x] = FindMatches(Map, KeyFrame,npoints)
X = [];
x = [];

if npoints > size(KeyFrame.ImagePoints,2)
    npoints = size(KeyFrame.ImagePoints,2);
end

for i = 1:size(KeyFrame.ImagePoints,2)
    input(1:3,i) = KeyFrame.ImagePoints(i).location;
    input(4,i) = KeyFrame.ImagePoints(i).id;
end

input = input(:,randperm(size(input,2)));

for i = 1:npoints
    for j = 1:size(Map.points,2)
        if (Map.points(j).id == input(4,i))
            X = [X Map.points(j).location];
            x = [x input(1:3,i)];
        end
    end
end


% --- Executes on button press in pushbutton_addkeyframe.
function pushbutton_addkeyframe_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_addkeyframe (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
AddKeyFrame(handles);

function AddKeyFrame(handles)
UpdateTick(handles);
PTAM = getappdata(handles.figure1,'ptam');
World = getappdata(handles.figure1,'world');
CurrKeyFrame = getappdata(handles.figure1,'currkeyframe');


PTAM.kfcount = PTAM.kfcount + 1;
PTAM.KeyFrames(PTAM.kfcount) = CurrKeyFrame;
DisplayKeyFrames(PTAM.KeyFrames, handles.viewtopdownest);

CurrKeyFrame.Camera = World.Camera;
World.KeyFrames(PTAM.kfcount) = CurrKeyFrame;


% Add some world points from this keyframe
KeyFrame1 = PTAM.KeyFrames(PTAM.kfcount);
KeyFrame2 =  FindClosestKeyFrame(PTAM.kfcount, PTAM.KeyFrames);
DisplayImage(KeyFrame1.ImagePoints, handles.viewkeyframe1);
DisplayImage(KeyFrame2.ImagePoints, handles.viewkeyframe2);

E1 = KeyFrame1.Camera.E(1:3,:);
E2 = KeyFrame2.Camera.E(1:3,:);
K = KeyFrame1.Camera.K;

outpoints = Reproject(KeyFrame1, KeyFrame2, E1, E2, K);
origsize = size(outpoints,2);
origmapsize = size(PTAM.Map.points,2);
PTAM.Map = AppendMap(PTAM.Map,outpoints);
newmapsize = size(PTAM.Map.points,2);
pointsadded = newmapsize-origmapsize;

display([int2str(pointsadded) ' points added out of ' int2str(origsize)]);
setappdata(handles.figure1,'ptam',PTAM);
setappdata(handles.figure1,'world',World);
UpdateTick(handles);


% --- Executes on button press in pushbutton_pathstep.
function pushbutton_pathstep_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_pathstep (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
theta = getappdata(handles.figure1,'theta');
World = getappdata(handles.figure1,'world');
theta = theta + 0.1;
ct = [18*cos(theta)-16 0 18*sin(theta)]';
World.Camera.camt = ct;
World.Camera.thetay = -theta;
World.Camera = RfromEuler(World.Camera);
setappdata(handles.figure1,'world',World);
setappdata(handles.figure1,'theta',theta);
EstimateCamera(handles);
UpdateTick(handles);


% --- Executes on button press in pushbutton_displaykf.
function pushbutton_displaykf_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_displaykf (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
PTAM = getappdata(handles.figure1,'ptam');
DisplayKeyFrames(PTAM.KeyFrames, handles.viewtopdownest);

function kfout = FindClosestKeyFrame(kfcount, KeyFrames)
mindist = 1000;
minkfcount = 0;


for i = 1:size(KeyFrames,2)
    if i ~= kfcount
        dist = norm(KeyFrames(i).Camera.E(1:3,4) - KeyFrames(kfcount).Camera.E(1:3,4));
        if dist < mindist
            mindist = dist;
            minkfcount = i;
        end
    end
end

kfout = KeyFrames(minkfcount);


function outMap = AppendMap(inMap,points)

outMap = inMap;

for i = 1:size(points,2)
    haspoint = false;
    for j = 1:size(inMap.points,2)
        if points(i).id == inMap.points(j).id
            haspoint = true;
        end
    end
    if haspoint == false
        outMap.points(size(outMap.points,2)+1) = points(i);
    end
    
    
end


% --- Executes on button press in pushbutton_bundle.
function pushbutton_bundle_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_bundle (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
PTAM = getappdata(handles.figure1,'ptam');
E1 = [];
for i = 1:PTAM.kfcount
    E1 = [E1; PTAM.KeyFrames(i).Camera.E];
end


PTAM = bundleadjust(PTAM);
setappdata(handles.figure1,'ptam',PTAM);
UpdateTick(handles);

E2 = [];
for i = 1:PTAM.kfcount
    E2 = [E2; PTAM.KeyFrames(i).Camera.E];
end

