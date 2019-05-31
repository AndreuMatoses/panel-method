 function varargout = GUI_Panels_method(varargin)

 %% GUI FOR THE PANEL METHOD IMPLEMENTATION
 % Authors:
 % Martínez Ibáñez, Josep Balbino
 % Matoses Gimenez, Andreu
 % Ortiz Moya, Álvaro
 % Ruiz López, Álvaro
 % Soriano Pascual, Fernando
 % 
 % UPV Aerodynamics
 % 31/05/2019
 
%% LICENCE
% MIT License Copyright (c) 2019 uerdnaGitHub
% 
% Permission is hereby granted, free of charge, to any person obtaining a
% copy of this software and associated documentation files (the
% "Software"), to deal in the Software without restriction, including
% without limitation the rights to use, copy, modify, merge, publish,
% distribute, sublicense, and/or sell copies of the Software, and to permit
% persons to whom the Software is furnished to do so, subject to the
% following conditions:
% 
% The above copyright notice and this permission notice shall be included
% in all copies or substantial portions of the Software.
% 
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
% OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
% MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN
% NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
% DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR
% OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE
% USE OR OTHER DEALINGS IN THE SOFTWARE.

 
%% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @GUI_Panels_method_OpeningFcn, ...
    'gui_OutputFcn',  @GUI_Panels_method_OutputFcn, ...
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


% --- Executes just before GUI_Panels_method is made visible.
function GUI_Panels_method_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to GUI_Panels_method (see VARARGIN)

% Choose default command line output for GUI_Panels_method
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes GUI_Panels_method wait for user response (see UIRESUME)
% uiwait(handles.figure1);

clear; clc;
% Define colors
global colors colorines 
colors.facemodel=[0.8 0.8 1.0]; colors.edgemodel='none'; colors.sectionplane='red';
colors.sectionpoints=[0.7,0,0.7]; colors.xysection='k'; colors.arrows='r';
colors.panelvertex=[0.7,0.2,0]; colors.panelmidpoints=[0.7,0.7,0];
colors.uppercp=[0.7,0.2,0]; colors.lowercp=[0.7,0.7,0];
colorines={[0 0 1]; [0 1 0]; [1 0 0]; [0 1 1]; [1 0 1]; [1 1 0]; [0 0 0]; [0.8 0.8 1.0]; [0.7,0,0.7]; [0.7,0.7,0]; [0.7,0.2,0]; [0 0 1]; [0 1 0]; [1 0 0]; [0 1 1]; [1 0 1]; [1 1 0]; [0 0 0]; [0.8 0.8 1.0]; [0.7,0,0.7]; [0.7,0.7,0]; [0.7,0.2,0]};
     
numberimage=1;
numberdata=1;
% Add path of used functions
addpath('functions\');


% --- Outputs from this function are returned to the command line.
function varargout = GUI_Panels_method_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in tglopen.
function tglopen_Callback(hObject, eventdata, handles)
% hObject    handle to tglopen (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global geometry section panels solutionplot colors solutionplot_fill solution_colorbar streamplot streamplot_fill contourplot contourplot_fill contourplot_colorbar 
% Reset all panels to zero
set(handles.txtsection,'String','Information about current section');
set(handles.tglsection,'Enable','off');
set(handles.lblsectionpoint,'Enable','off');
set(handles.lblsectionvector,'Enable','off');
set(handles.tgldiscretize,'Enable','off');
set(handles.lblnpanels,'Enable','off');
set(handles.lblslpsens,'Enable','off');
set(handles.poptypesolver,'Enable','off');
set(handles.lblangleincidencestart,'Enable','off');
set(handles.lblangleincidencestep,'Enable','off');
set(handles.lblangleincidenceend,'Enable','off');
set(handles.tglsolve,'Enable','off');
set(handles.lblvelocity,'Enable','off');
set(handles.tglinvert,'Enable','off');
set(handles.lblrotate,'Enable','off');
set(handles.lblxtrans,'Enable','off');
set(handles.lblytrans,'Enable','off');
% File explorer to select geometry
[geometry.name, geometry.path]=uigetfile('.stl','Select a geometry');
% Check that an geometry in .stl format has been selected
if strfind(geometry.name,'.stl')
    % Indicate the name of the geometry file in the static text: txtopen
    set(handles.txtopen,'String',geometry.name);
    % Add the path of the geometry and read it
    addpath(geometry.path);
    geometry.model=stlread(geometry.name);
    % Delete previous plot and plot the new geometry, lights, material...
    try
        delete(geometry.plot.model)
        delete(geometry.plot.light)
    catch
    end
    try
        delete(section.plot.surf)
        delete(section.plot.intersec)
    catch
    end
    try
        delete(section.plot.xyintersec)
        delete(section.plot.arrows)
        delete(panels.plot)
        delete(panels.plot_text)
        delete(panels.plot_midpoints)
    catch
    end
    try
        delete(solutionplot)
        delete(solutionplot_fit)
    end
    try
        delete(streamplot)
        delete(streamplot_fill)
    end
    try
        delete(contourplot)
        delete(contourplot_fill)
        delete(contourplot_colorbar)
    end
    geometry.plot.model=patch(handles.axplot,geometry.model,'FaceColor',colors.facemodel,...
        'EdgeColor',colors.edgemodel,'FaceLighting','gouraud',...
        'AmbientStrength',0.15,'FaceAlpha',0.5); % Represent the model
    hold(handles.axplot,'on'); % Hold the graph to add other geometries
    geometry.plot.light=light(handles.axplot); % Add a light
    xlabel(handles.axplot,'x [m]'); % Label of each axes
    ylabel(handles.axplot,'y [m]');
    zlabel(handles.axplot,'z [m]');
    title('')
    material('metal'); % Change the material of the model
    axis('image'); % Select the scale of the axes
    view([45 35]); % Indicate the perspective angle
    legend('Geometry','Location','northeast'); % Add the legend
    % Activate the selection of the section plane
    set(handles.lblsectionpoint,'Enable','on');
    set(handles.lblsectionvector,'Enable','on');
    % Read the position of the section plane and plot it
    section_evaluation(hObject, eventdata, handles);
% If the selected file is not a valid .stl file indicate it
else
    % Delete, if any, geometry and section
    try
        delete(geometry.plot.model)
        delete(geometry.plot.light)
    catch
    end
    try
        delete(section.plot.surf)
        delete(section.plot.intersec)
    catch
    end
    try
        delete(section.plot.xyintersec)
        delete(section.plot.arrows)
        delete(panels.plot)
        delete(panels.plot_text)
        delete(panels.plot_midpoints)
    catch
    end
    try
        delete(solutionplot)
        delete(solutionplot_fit)
    end
    try
        delete(streamplot)
        delete(streamplot_fill)
    end
    try
        delete(contourplot)
        delete(contourplot_fill)
        delete(contourplot_colorbar)
    end
    % Change texts and status of buttons
    set(handles.txtopen,'String','Invalid file');
    set(handles.txtsection,'String','Information about current section');
    set(handles.tglsection,'Enable','off');
    set(handles.lblsectionpoint,'Enable','off');
    set(handles.lblsectionvector,'Enable','off');
    set(handles.tgldiscretize,'Enable','off');
    set(handles.lblnpanels,'Enable','off');
    set(handles.lblslpsens,'Enable','off');
    set(handles.poptypesolver,'Enable','off');
    set(handles.lblangleincidencestart,'Enable','off');
    set(handles.lblangleincidencestep,'Enable','off');
    set(handles.lblangleincidenceend,'Enable','off');
    set(handles.tglsolve,'Enable','off');
    set(handles.lblvelocity,'Enable','off');
    set(handles.tglinvert,'Enable','off');
    set(handles.lblrotate,'Enable','off');
    set(handles.lblxtrans,'Enable','off');
    set(handles.lblytrans,'Enable','off');
end


function lblsectionpoint_Callback(hObject, eventdata, handles)
% hObject    handle to lblsectionpoint (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of lblsectionpoint as text
%        str2double(get(hObject,'String')) returns contents of lblsectionpoint as a double
section_evaluation(hObject, eventdata, handles);


% --- Executes during object creation, after setting all properties.
function lblsectionpoint_CreateFcn(hObject, eventdata, handles)
% hObject    handle to lblsectionpoint (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function lblsectionvector_Callback(hObject, eventdata, handles)
% hObject    handle to lblsectionvector (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of lblsectionvector as text
%        str2double(get(hObject,'String')) returns contents of lblsectionvector as a double
section_evaluation(hObject, eventdata, handles);


% --- Executes during object creation, after setting all properties.
function lblsectionvector_CreateFcn(hObject, eventdata, handles)
% hObject    handle to lblsectionvector (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Function to calculate the section and represent it in the plot
function section_evaluation(hObject, eventdata, handles)
% This function is evaluated either when a geometry is selected to check if
% default section is valid or when the user manually changes the position
% or orientation of the section.
global geometry section colors 
% Delete, if any, previous sections
try
    delete(section.plot.surf)
    delete(section.plot.intersec)
catch
end
% Read the origin where the section will be placed
section.origin=str2num(get(handles.lblsectionpoint,'String'));
% Check that this points has 3 coordinates, if not, set an empty vector
if length(section.origin)~=3
    section.origin=[];
end
% Read the orientation of the section plane
section.normal=get(handles.lblsectionvector,'String');
% Check that it is a valid orientation, if not, set an empty vector
if strlength(section.normal)==1
    if strfind(section.normal,'x')
        section.normal=[1,0,0];
    elseif strfind(section.normal,'y')
        section.normal=[0,1,0];
    elseif strfind(section.normal,'z')
        section.normal=[0,0,1];
    else
        section.normal=[];
    end
else
    section.normal=[];
end
% If the origin and the orientation of the section plane, represent it
if (not(isempty(section.origin))&&not(isempty(section.normal)))
    % Find geometry points enclosed by the section
    section.valid_indexes=geometry.model.vertices(:,section.normal==1)==section.origin(section.normal==1);
    section.points=geometry.model.vertices(section.valid_indexes,:);
    % Delete repeated points of this intersection
    section.points=unique(section.points,'rows');
    % Sort these points in a closed pattern
    section=angular_sorting(section);
    % If the intersection is not empty
    if not(isempty(section.points))
        % Find two orthonormal vectors which are orthogonal to v
        section.ortovectors = null(section.normal);
        section.max_points=max(section.points);
        section.min_points=min(section.points);
        section.max_points=section.max_points(section.normal==0);
        section.min_points=section.min_points(section.normal==0);
        % Provide a gridwork (you choose the size)
        [section.P,section.Q] = meshgrid(-[section.max_points(1),section.min_points(1)],[section.max_points(2),section.min_points(2)]);
        % Compute the corresponding cartesian coordinates
        % using the two vectors in section.ortovectors
        section.X = section.origin(1)+section.ortovectors(1,1)*section.P+section.ortovectors(1,2)*section.Q;
        section.Y = section.origin(2)+section.ortovectors(2,1)*section.P+section.ortovectors(2,2)*section.Q;
        section.Z = section.origin(3)+section.ortovectors(3,1)*section.P+section.ortovectors(3,2)*section.Q;
        section.plot.surf=surf(section.X,section.Y,section.Z,'FaceAlpha',0.5,'FaceColor',colors.sectionplane);
        % Plot intersection points
        section.plot.intersec=plot3(handles.axplot,section.points(:,1),section.points(:,2),section.points(:,3),'*','MarkerEdgeColor',colors.sectionpoints);
        legend('Geometry','Section','Intersection points','Location','northeast'); % Add the legend
        % Indicate that the section is valid and activate tglsection
        set(handles.txtsection,'String','Valid section');
        set(handles.tglsection,'Enable','on');
    % If the intersection is not valid, indicate it and desactivate tglsection
    else
        set(handles.txtsection,'String','Invalid section');
        set(handles.tglsection,'Enable','off');
    end
% If the intersection is not valid, indicate it and desactivate tglsection
else
    set(handles.txtsection,'String','Invalid section');
    set(handles.tglsection,'Enable','off');
end


% --- Executes on button press in tglsection.
function tglsection_Callback(hObject, eventdata, handles)
% hObject    handle to tglsection (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of tglsection
global section geometry colors 
% Prepare the points in two dimensions
section.xy_points=section.points(:,not(section.normal));
% Delete the plot to represent only the section
try
    delete(geometry.plot.model)
    delete(geometry.plot.light)
catch
end
try
    delete(section.plot.surf)
    delete(section.plot.intersec)
catch
end
% Represent in the graph the two dimensional intersection
section.plot.xyintersec=plot(handles.axplot,[section.xy_points(:,1);section.xy_points(1,1)],[section.xy_points(:,2);section.xy_points(1,2)],colors.xysection);
view([0,90]); % Set perspective of the intersection, top view
% Enable the discretization panel and section options
set(handles.tgldiscretize,'Enable','on');
set(handles.lblnpanels,'Enable','on');
set(handles.lblslpsens,'Enable','on');
set(handles.tglinvert,'Enable','on');
set(handles.lblrotate,'Enable','on');
set(handles.lblxtrans,'Enable','on');
set(handles.lblytrans,'Enable','on');
% Desactivate the section panel and indicate it
set(handles.tglsection,'Value',0);
set(handles.tglsection,'Enable','off');
set(handles.lblsectionpoint,'Enable','off');
set(handles.lblsectionvector,'Enable','off');
set(handles.txtsection,'String','Section correctly selected')
% Plot the free stream
section.arrowY=linspace(min(section.xy_points(:,2)),max(section.xy_points(:,2)),5);
section.arrowX=ones(size(section.arrowY))*(min(section.xy_points(:,1))-(max(section.xy_points(:,1))-min(section.xy_points(:,1)))/5);
section.arrowDX=ones(size(section.arrowY))*(max(section.xy_points(:,1))-min(section.xy_points(:,1)))/7;
section.arrowDY=ones(size(section.arrowY))*0;
section.plot.arrows=quiver(section.arrowX,section.arrowY,section.arrowDX,section.arrowDY,0,colors.arrows);
legend('Section','Freestream','Location','southwest','Location','southwest'); % Add the legend


function lblnpanels_Callback(hObject, eventdata, handles)
% hObject    handle to lblnpanels (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of lblnpanels as text
%        str2double(get(hObject,'String')) returns contents of lblnpanels as a double


% --- Executes during object creation, after setting all properties.
function lblnpanels_CreateFcn(hObject, eventdata, handles)
% hObject    handle to lblnpanels (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function lblslpsens_Callback(hObject, eventdata, handles)
% hObject    handle to lblslpsens (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of lblslpsens as text
%        str2double(get(hObject,'String')) returns contents of lblslpsens as a double


% --- Executes during object creation, after setting all properties.
function lblslpsens_CreateFcn(hObject, eventdata, handles)
% hObject    handle to lblslpsens (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in tgldiscretize.
function tgldiscretize_Callback(hObject, eventdata, handles)
% hObject    handle to tgldiscretize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global section panels solutionplot colors solutionplot_fit streamplot streamplot_fill contourplot contourplot_fill contourplot_colorbar 
% Desactivate section options 
set(handles.tglinvert,'Enable','off');
set(handles.lblrotate,'Enable','off');
set(handles.lblxtrans,'Enable','off');
set(handles.lblytrans,'Enable','off');
% Delete previous sections, arrows and discretizations, if any
try 
    delete(section.plot.xyintersec)
    delete(section.plot.arrows)
    delete(panels.plot)
    delete(panels.plot_text)
    delete(panels.plot_midpoints)
catch
end
%Delete previous solutions if any
try 
    delete(solutionplot)
    delete(solutionplot_fit)
catch
end
try
    delete(streamplot)
    delete(streamplot_fill)
end
try
    delete(contourplot)
    delete(contourplot_fill)
    delete(contourplot_colorbar)
end
% Label of the axes
xlabel(handles.axplot,'x [m]')
ylabel(handles.axplot,'y [m]')
title('')
% Replot the section and the free stream
section.plot.xyintersec=plot(handles.axplot,[section.xy_points(:,1);section.xy_points(1,1)],[section.xy_points(:,2);section.xy_points(1,2)],colors.xysection);
section.plot.arrows=quiver(section.arrowX,section.arrowY,section.arrowDX,section.arrowDY,0,colors.arrows);
panels.npanels=str2double(get(handles.lblnpanels,'String')); % Get number of panels desired by the user
panels.slpsens=str2double(get(handles.lblslpsens,'String')); % Get slope sensitivity desired by the user
% According to selected values calculate the panels
panels=panel_generator(section.xy_points,panels.npanels,panels.slpsens);
% Plot the panels
panels.plot=plot(handles.axplot,[panels.vertex(:,1); panels.vertex(1,1)],[panels.vertex(:,2); panels.vertex(1,2)],'-o','Color',colors.panelvertex,'MarkerEdgeColor',colors.panelvertex);
panels.plot_midpoints=plot(handles.axplot,panels.mid_points(:,1),panels.mid_points(:,2),'*','MarkerEdgeColor',colors.panelmidpoints);
panels.plot_text=text(panels.mid_points(:,1),panels.mid_points(:,2),num2str([1:length(panels.mid_points(:,2))]'));
legend('Section','Freestream','Panels','Midpoints','Location','southwest'); % Add the legend
axis('image')
% Enable the solution button
set(handles.poptypesolver,'Enable','on');
set(handles.tglsolve,'Enable','on');
set(handles.lblvelocity,'Enable','on');
if (get(handles.poptypesolver,'Value')==1)
    set(handles.lblangleincidencestart,'Enable','off');
    set(handles.lblangleincidencestep,'Enable','off');
    set(handles.lblangleincidenceend,'Enable','off');
elseif (get(handles.poptypesolver,'Value')==2)
    set(handles.lblangleincidencestart,'Enable','on');
    set(handles.lblangleincidencestep,'Enable','on');
    set(handles.lblangleincidenceend,'Enable','on');
end


function lblvelocity_Callback(hObject, eventdata, handles)
% hObject    handle to lblvelocity (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of lblvelocity as text
%        str2double(get(hObject,'String')) returns contents of lblvelocity as a double
set(handles.tglsolve,'Enable','on')

% --- Executes during object creation, after setting all properties.
function lblvelocity_CreateFcn(hObject, eventdata, handles)
% hObject    handle to lblvelocity (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in tglsolve.
function tglsolve_Callback(hObject, eventdata, handles)
% hObject    handle to tglsolve (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of tglsolve
global panels solution section solutionplot flow solutionplot_fit streamplot streamplot_fill  contourplot contourplot_fill contourplot_colorbar 
% Get desired velocity of the free stream by the user
flow.velocity_infinite=str2double(get(handles.lblvelocity,'String'));
% Solve the proble considering the type of the solver
if (get(handles.poptypesolver,'Value')==1)
    flow.type=1;
    flow.alpha=0; %alpha in deg 
    solucion(1)=panel_solveNL(panels,flow.velocity_infinite,flow.alpha);
    solution=solucion;
elseif (get(handles.poptypesolver,'Value')==2)
    flow.type=2;
    flow.alpha=linspace(str2double(get(handles.lblangleincidencestart,'String')),...
        str2double(get(handles.lblangleincidenceend,'String')),str2double(get(handles.lblangleincidencestep,'String')));
    for i=1:length(flow.alpha)
        solucion(i)=panel_solverL_iter(panels,flow.velocity_infinite,flow.alpha(i));
    end
    solution=solucion;
end
% Delete previous solutions or discretizations
try
    delete(section.plot.xyintersec)
    delete(section.plot.arrows)
    delete(panels.plot)
    delete(panels.plot_text)
    delete(panels.plot_midpoints)
catch
end
try
    delete(solutionplot)
    delete(solutionplot_fit)
end
try
    delete(streamplot)
    delete(streamplot_fill)
end
try
    delete(contourplot)
    delete(contourplot_fill)
    delete(contourplot_colorbar)
end
axis('normal')
% Indicate which angles of incidence are calculated
flow.txtalpha=[];
for i=1:length(flow.alpha)
    flow.txtalpha{i,1}=num2str(flow.alpha(i));
end
% Enable post processing panel
set(handles.popsolutions,'String',flow.txtalpha);
set(handles.chkcpplot,'Enable','on');
set(handles.chkclplot,'Enable','on');
set(handles.chkstreamlines,'Enable','on');
set(handles.chkpressure,'Enable','on');
set(handles.lblresolution,'Enable','on');
set(handles.lblsizex,'Enable','on');
set(handles.lblsizey,'Enable','on');
set(handles.lblnstreamlines,'Enable','on');
set(handles.popxyabs,'Enable','on');
set(handles.chkvelocitycontour,'Enable','on');
% Indicate that the cp plot is the one selected, and indicate the first angle of incidence
set(handles.chkcpplot,'Value',1);
set(handles.chkclplot,'Value',0);
flow.alphaplot=zeros(size(flow.alpha));
flow.alphaplot(1)=1;
set(handles.txtcpplot,'String', ['alphas: ' num2str(flow.alpha(1))])
% Plot the solution
solutionplot=plot(panels.upperx,solution(1).Cpupper,'-^',panels.lowerx,solution(1).Cplower,'-o','Color',[0,0,1]);
legend(['Upper Alpha= ' num2str(flow.alpha(1))],['Lower Alpha= ' num2str(flow.alpha(1))],'Location','southwest'); % Add the legend
ylabel('Pressure coefficient, -Cp');
xlabel('x/c');
title('')
% Disable the solution button
set(handles.tglsolve,'Enable','off');


function lblrotate_Callback(hObject, eventdata, handles)
% hObject    handle to lblrotate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of lblrotate as text
%        str2double(get(hObject,'String')) returns contents of lblrotate as a double
global section colors 

% Delete previous section and free stream plots
try 
    delete(section.plot.xyintersec)
    delete(section.plot.arrows)
catch
end
% Calculate the new positions of the points and plot them 
section.rot_alfa=str2double(get(hObject,'String'));
section.xy_points=section.xy_points*[cosd(section.rot_alfa),-sind(section.rot_alfa);sind(section.rot_alfa),cosd(section.rot_alfa)];
section.plot.xyintersec=plot(handles.axplot,[section.xy_points(:,1);section.xy_points(1,1)],[section.xy_points(:,2);section.xy_points(1,2)],colors.xysection);
% Plot again the free stream
section.arrowY=linspace(min(section.xy_points(:,2)),max(section.xy_points(:,2)),5);
section.arrowX=ones(size(section.arrowY))*(min(section.xy_points(:,1))-(max(section.xy_points(:,1))-min(section.xy_points(:,1)))/5);
section.arrowDX=ones(size(section.arrowY))*(max(section.xy_points(:,1))-min(section.xy_points(:,1)))/7;
section.arrowDY=ones(size(section.arrowY))*0;
section.plot.arrows=quiver(section.arrowX,section.arrowY,section.arrowDX,section.arrowDY,0,colors.arrows);
legend('Section','Freestream','Location','southwest'); % Add the legend


% --- Executes during object creation, after setting all properties.
function lblrotate_CreateFcn(hObject, eventdata, handles)
% hObject    handle to lblrotate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function lblxtrans_Callback(hObject, eventdata, handles)
% hObject    handle to lblxtrans (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of lblxtrans as text
%        str2double(get(hObject,'String')) returns contents of lblxtrans as a double
global section colors 

% Delete previous section and free stream plots
try 
    delete(section.plot.xyintersec)
    delete(section.plot.arrows)
catch
end
% Calculate the new positions of the points and plot them 
section.x_trans=str2double(get(hObject,'String'));
section.xy_points=[section.xy_points(:,1)+section.x_trans,section.xy_points(:,2)];
section.plot.xyintersec=plot(handles.axplot,[section.xy_points(:,1);section.xy_points(1,1)],[section.xy_points(:,2);section.xy_points(1,2)],colors.xysection);
% Plot again the free stream
section.arrowY=linspace(min(section.xy_points(:,2)),max(section.xy_points(:,2)),5);
section.arrowX=ones(size(section.arrowY))*(min(section.xy_points(:,1))-(max(section.xy_points(:,1))-min(section.xy_points(:,1)))/5);
section.arrowDX=ones(size(section.arrowY))*(max(section.xy_points(:,1))-min(section.xy_points(:,1)))/7;
section.arrowDY=ones(size(section.arrowY))*0;
section.plot.arrows=quiver(section.arrowX,section.arrowY,section.arrowDX,section.arrowDY,0,colors.arrows);
legend('Section','Freestream','Location','southwest'); % Add the legend


% --- Executes during object creation, after setting all properties.
function lblxtrans_CreateFcn(hObject, eventdata, handles)
% hObject    handle to lblxtrans (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function lblytrans_Callback(hObject, eventdata, handles)
% hObject    handle to lblytrans (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of lblytrans as text
%        str2double(get(hObject,'String')) returns contents of lblytrans as a double
global section colors 

% Delete previous section and free stream plots
try 
    delete(section.plot.xyintersec)
    delete(section.plot.arrows)
catch
end
% Calculate the new positions of the points and plot them 
section.y_trans=str2double(get(hObject,'String'));
section.xy_points=[section.xy_points(:,1),section.xy_points(:,2)+section.y_trans];
section.plot.xyintersec=plot(handles.axplot,[section.xy_points(:,1);section.xy_points(1,1)],[section.xy_points(:,2);section.xy_points(1,2)],colors.xysection);
% Plot again the free stream
section.arrowY=linspace(min(section.xy_points(:,2)),max(section.xy_points(:,2)),5);
section.arrowX=ones(size(section.arrowY))*(min(section.xy_points(:,1))-(max(section.xy_points(:,1))-min(section.xy_points(:,1)))/5);
section.arrowDX=ones(size(section.arrowY))*(max(section.xy_points(:,1))-min(section.xy_points(:,1)))/7;
section.arrowDY=ones(size(section.arrowY))*0;
section.plot.arrows=quiver(section.arrowX,section.arrowY,section.arrowDX,section.arrowDY,0,colors.arrows);
legend('Section','Freestream','Location','southwest'); % Add the legend


% --- Executes during object creation, after setting all properties.
function lblytrans_CreateFcn(hObject, eventdata, handles)
% hObject    handle to lblytrans (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in tglinvert.
function tglinvert_Callback(hObject, eventdata, handles)
% hObject    handle to tglinvert (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global section colors 

% Delete previous section and free stream plots
try 
    delete(section.plot.xyintersec)
    delete(section.plot.arrows)
catch
end
% Calculate the new positions of the points and plot them 
section.xy_points=[-section.xy_points(:,1),section.xy_points(:,2)];
section.plot.xyintersec=plot(handles.axplot,[section.xy_points(:,1);section.xy_points(1,1)],[section.xy_points(:,2);section.xy_points(1,2)],colors.xysection);
% Plot again the free stream
section.arrowY=linspace(min(section.xy_points(:,2)),max(section.xy_points(:,2)),5);
section.arrowX=ones(size(section.arrowY))*(min(section.xy_points(:,1))-(max(section.xy_points(:,1))-min(section.xy_points(:,1)))/5);
section.arrowDX=ones(size(section.arrowY))*(max(section.xy_points(:,1))-min(section.xy_points(:,1)))/7;
section.arrowDY=ones(size(section.arrowY))*0;
section.plot.arrows=quiver(section.arrowX,section.arrowY,section.arrowDX,section.arrowDY,0,colors.arrows);
legend('Section','Freestream','Location','southwest'); % Add the legend


% --- Executes on button press in tglpan.
function tglpan_Callback(hObject, eventdata, handles)
% hObject    handle to tglpan (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
rotate3d(handles.axplot,'off'); % Disable manual rotation of the graph
zoom(handles.axplot,'off'); % Disable zoom in the plot
pan(handles.axplot,'on'); % Enable manual pan of the graph

% --- Executes on button press in tglrotate.
function tglrotate_Callback(hObject, eventdata, handles)
% hObject    handle to tglrotate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
rotate3d(handles.axplot,'on'); % Enable manual rotation of the graph
zoom(handles.axplot,'off'); % Disable zoom in the plot
pan(handles.axplot,'off'); % Disable manual pan of the graph

% --- Executes on button press in tglzoom.
function tglzoom_Callback(hObject, eventdata, handles)
% hObject    handle to tglzoom (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
rotate3d(handles.axplot,'off'); % Disable manual rotation of the graph
zoom(handles.axplot,'on'); % Enable zoom in the plot
pan(handles.axplot,'off'); % Disable manual pan of the graph

% --- Executes on button press in tglmouse.
function tglmouse_Callback(hObject, eventdata, handles)
% hObject    handle to tglmouse (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
rotate3d(handles.axplot,'off'); % Disable manual rotation of the graph
zoom(handles.axplot,'off'); % Disable zoom in the plot
pan(handles.axplot,'off'); % Disable manual pan of the graph

% --- Executes on selection change in poptypesolver.
function poptypesolver_Callback(hObject, eventdata, handles)
% hObject    handle to poptypesolver (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns poptypesolver contents as cell array
%        contents{get(hObject,'Value')} returns selected item from
%        poptypesolver
global solutionplot solutionplot_fit streamplot streamplot_fill contourplot contourplot_fill contourplot_colorbar 
if (get(handles.poptypesolver,'Value')==1)
    set(handles.lblangleincidencestart,'Enable','off');
    set(handles.lblangleincidencestep,'Enable','off');
    set(handles.lblangleincidenceend,'Enable','off');
    set(handles.tglsolve,'Enable','on')
elseif (get(handles.poptypesolver,'Value')==2)
    set(handles.lblangleincidencestart,'Enable','on');
    set(handles.lblangleincidencestep,'Enable','on');
    set(handles.lblangleincidenceend,'Enable','on');
    set(handles.tglsolve,'Enable','on')
end
try
    delete(solutionplot)
    delete(solutionplot_fit)
catch
end
try
    delete(streamplot)
    delete(streamplot_fill)
end
try
    delete(contourplot)
    delete(contourplot_fill)
    delete(contourplot_colorbar)
end
    
% --- Executes during object creation, after setting all properties.
function poptypesolver_CreateFcn(hObject, eventdata, handles)
% hObject    handle to poptypesolver (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function lblangleincidencestart_Callback(hObject, eventdata, handles)
% hObject    handle to lblangleincidencestart (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of lblangleincidencestart as text
%        str2double(get(hObject,'String')) returns contents of lblangleincidencestart as a double
set(handles.tglsolve,'Enable','on')

% --- Executes during object creation, after setting all properties.
function lblangleincidencestart_CreateFcn(hObject, eventdata, handles)
% hObject    handle to lblangleincidencestart (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function lblangleincidencestep_Callback(hObject, eventdata, handles)
% hObject    handle to lblangleincidencestep (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of lblangleincidencestep as text
%        str2double(get(hObject,'String')) returns contents of lblangleincidencestep as a double
set(handles.tglsolve,'Enable','on')

% --- Executes during object creation, after setting all properties.
function lblangleincidencestep_CreateFcn(hObject, eventdata, handles)
% hObject    handle to lblangleincidencestep (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function lblangleincidenceend_Callback(hObject, eventdata, handles)
% hObject    handle to lblangleincidenceend (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of lblangleincidenceend as text
%        str2double(get(hObject,'String')) returns contents of lblangleincidenceend as a double
set(handles.tglsolve,'Enable','on')

% --- Executes during object creation, after setting all properties.
function lblangleincidenceend_CreateFcn(hObject, eventdata, handles)
% hObject    handle to lblangleincidenceend (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popsolutions.
function popsolutions_Callback(hObject, eventdata, handles)
% hObject    handle to popsolutions (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popsolutions contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popsolutions


% --- Executes during object creation, after setting all properties.
function popsolutions_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popsolutions (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in chkcpplot.
function chkcpplot_Callback(hObject, eventdata, handles)
% hObject    handle to chkcpplot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of chkcpplot
global solutionplot flow colorines panels solution solutionplot_fit streamplot streamplot_fill contourplot contourplot_fill contourplot_colorbar 
flow.alphaplotval=flow.alpha(flow.alphaplot==1);
flow.numberalphas=length(flow.alphaplotval);
if (get(hObject,'Value'))
    try
        delete(solutionplot)
        delete(solutionplot_fit)
    end
    try
        delete(streamplot)
        delete(streamplot_fill)
    end
    try
        delete(contourplot)
        delete(contourplot_fill)
        delete(contourplot_colorbar)
    end
    set(handles.chkclplot,'Value',0);
    % Colors of the plot
    indexes=[1:flow.numberalphas,12:(12+flow.numberalphas-1)];
    for i=1:length(indexes)
        colores{i,1}=colorines{indexes(i),1};
    end
    % Plot the solutions
    solutionplot=plot(panels.upperx,[solution(flow.alphaplot==1).Cpupper],'-^',panels.lowerx,[solution(flow.alphaplot==1).Cplower],'-o');
    set(solutionplot, {'color'}, colores);
    ylabel('Pressure coefficient, -Cp');
    xlabel('x/c');
    axis('square');
    title('')
    leg={};
    for i=1:flow.numberalphas
        leg{i}=['Upper Alpha= ' num2str(flow.alphaplotval(i))];
        leg{i+flow.numberalphas}=['Lower Alpha= ' num2str(flow.alphaplotval(i))];
    end
    legend(handles.axplot,leg);
else
    % Delete previous solution plot
    try
        delete(solutionplot)
        delete(solutionplot_fit)
    end
    try
        delete(streamplot)
        delete(streamplot_fill)
    end
    try
        delete(contourplot)
        delete(contourplot_fill)
        delete(contourplot_colorbar)
    end

end


% --- Executes on button press in tgladdplot.
function tgladdplot_Callback(hObject, eventdata, handles)
% hObject    handle to tgladdplot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global flow solutionplot panels solution colorines solutionplot_fit streamplot streamplot_fill contourplot contourplot_fill contourplot_colorbar 
% Get the selected alpha case to represent
lectura=get(handles.popsolutions,'Value');
flow.alphaplot(lectura)=not(flow.alphaplot(lectura));
flow.alphaplotval=flow.alpha(flow.alphaplot==1);
flow.numberalphas=length(flow.alphaplotval);
% Delete previous solution plot
try
    delete(solutionplot)
    delete(solutionplot_fit)
end
try
    delete(streamplot)
    delete(streamplot_fill)
end
try
    delete(contourplot)
    delete(contourplot_fill)
    delete(contourplot_colorbar)
end

% Show the alphas in the plot
if isempty(flow.alphaplotval)
    set(handles.txtcpplot,'String','No case selected')
else
    set(handles.txtcpplot,'String', ['alphas: ' mat2str(round(flow.alphaplotval,1))])
    if get(handles.chkcpplot,'Value')==1
    % Colors of the plot
    indexes=[1:flow.numberalphas,12:(12+flow.numberalphas-1)];
    for i=1:length(indexes)
        colores{i,1}=colorines{indexes(i),1};
    end
    % Plot the solutions
    solutionplot=plot(panels.upperx,[solution(flow.alphaplot==1).Cpupper],'-^',panels.lowerx,[solution(flow.alphaplot==1).Cplower],'-o');
    set(solutionplot, {'color'}, colores);
    axis('square'); % Select the scale of the axes
    ylabel('Pressure coefficient, -Cp');
    xlabel('x/c');
   title('')
    leg={};
    for i=1:flow.numberalphas
        leg{i}=['Upper Alpha= ' num2str(flow.alphaplotval(i))];
        leg{i+flow.numberalphas}=['Lower Alpha= ' num2str(flow.alphaplotval(i))];
    end
    legend(handles.axplot,leg);
    elseif get(handles.chkclplot,'Value')==1
        solutionplot=plot(handles.axplot,flow.alphaplotval,[solution(flow.alphaplot==1).cl_kutta],'-o');
        coeffs = polyfit(flow.alphaplotval, [solution(flow.alphaplot==1).cl_kutta], 1);
        xFitting = flow.alphaplotval;
        yFitted = polyval(coeffs, xFitting);      
        solutionplot_fit=plot(xFitting, yFitted, 'r', 'LineWidth', 2);
        xlabel('alpha [º]');
        ylabel('cl Kutta-Joukowski');    
        legend('Calculated values', ['Least square fit: Slope ' num2str(coeffs(1)) ', Intercept: ' num2str(coeffs(2))],'Location','southeast');
        grid on
    elseif get(handles.chkcmplot,'Value')==1
    end
end



% --- Executes on button press in chkclplot.
function chkclplot_Callback(hObject, eventdata, handles)
% hObject    handle to chkclplot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of chkclplot
global solutionplot flow colorines panels solution solutionplot_fit streamplot streamplot_fill contourplot contourplot_fill contourplot_colorbar 
flow.alphaplotval=flow.alpha(flow.alphaplot==1);
flow.numberalphas=length(flow.alphaplotval);
if (get(hObject,'Value'))
    try
        delete(solutionplot)
        delete(solutionplot_fit)
    end
    try
        delete(streamplot)
        delete(streamplot_fill)
    end
    try
        delete(contourplot)
        delete(contourplot_fill)
        delete(contourplot_colorbar)
    end
    set(handles.chkcpplot,'Value',0);
    solutionplot=plot(handles.axplot,flow.alphaplotval,[solution(flow.alphaplot==1).cl_kutta],'-o');
    coeffs = polyfit(flow.alphaplotval, [solution(flow.alphaplot==1).cl_kutta], 1);
    xFitting = flow.alphaplotval;
    yFitted = polyval(coeffs, xFitting);
    solutionplot_fit=plot(xFitting, yFitted, 'r', 'LineWidth', 2);
    xlabel('alpha [º]');
    ylabel('cl Kutta-Joukowski');
    legend('Calculated values', ['Least square fit: Slope ' num2str(coeffs(1)) ', Intercept: ' num2str(coeffs(2))],'Location','southeast');
    grid on
else
    % Delete previous solution plot
    try
        delete(solutionplot)
        delete(solutionplot_fit)
    end
    try
        delete(streamplot)
        delete(streamplot_fill)
    end
    try
        delete(contourplot)
        delete(contourplot_fill)
        delete(contourplot_colorbar)
    end

end

% --- Executes on button press in chkcmplot.
function chkcmplot_Callback(hObject, eventdata, handles)
% hObject    handle to chkcmplot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of chkcmplot
global solutionplot solutionplot_fit streamplot streamplot_fill contourplot contourplot_fill contourplot_colorbar 
if (get(hObject,'Value'))
    set(handles.chkcpplot,'Value',0);
    set(handles.chkclplot,'Value',0);
else
    % Delete previous solution plot
    try
        delete(solutionplot)
        delete(solutionplot_fit)
    end
    try
        delete(streamplot)
        delete(streamplot_fill)
    end
    try
        delete(contourplot)
        delete(contourplot_fill)
        delete(contourplot_colorbar)
    end
end


% --- Executes on button press in tgladdallplot.
function tgladdallplot_Callback(hObject, eventdata, handles)
% hObject    handle to tgladdallplot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global flow solutionplot panels solution colorines streamplot streamplot_fill solutionplot_fit contourplot contourplot_fill contourplot_colorbar 
% See if represent or delete all alphas 
flow.alphaplotval=flow.alpha(flow.alphaplot==1);
flow.numberalphas=length(flow.alphaplotval);
if flow.numberalphas<length(flow.alpha)
    flow.alphaplot=ones([1,length(flow.alpha)]);
else
    flow.alphaplot=zeros([1,length(flow.alpha)]);
end
flow.alphaplotval=flow.alpha(flow.alphaplot==1);
flow.numberalphas=length(flow.alphaplotval);
% Delete previous solution plot
try
    delete(solutionplot)
    delete(solutionplot_fit)
end
try
    delete(streamplot)
    delete(streamplot_fill)
end
try
    delete(contourplot)
    delete(contourplot_fill)
    delete(contourplot_colorbar)
end

% Show the alphas in the plot
if isempty(flow.alphaplotval)
    set(handles.txtcpplot,'String','No case selected')
else
    set(handles.txtcpplot,'String', ['alphas: ' mat2str(round(flow.alphaplotval,1))])
    if get(handles.chkcpplot,'Value')==1
        % Colors of the plot
        indexes=[1:flow.numberalphas,12:(12+flow.numberalphas-1)];
        for i=1:length(indexes)
            colores{i,1}=colorines{indexes(i),1};
        end
        % Plot the solutions
        solutionplot=plot(panels.upperx,[solution(flow.alphaplot==1).Cpupper],'-^',panels.lowerx,[solution(flow.alphaplot==1).Cplower],'-o');
        set(solutionplot, {'color'}, colores);
        axis('square'); % Select the scale of the axes
        ylabel('Pressure coefficient, -Cp');
        xlabel('x/c');
        title('')
        leg={};
        for i=1:flow.numberalphas
            leg{i}=['Upper Alpha= ' num2str(flow.alphaplotval(i))];
            leg{i+flow.numberalphas}=['Lower Alpha= ' num2str(flow.alphaplotval(i))];
        end
        legend(handles.axplot,leg);
    elseif get(handles.chkclplot,'Value')==1
        solutionplot=plot(handles.axplot,flow.alphaplotval,[solution(flow.alphaplot==1).cl_kutta],'-o');
        coeffs = polyfit(flow.alphaplotval, [solution(flow.alphaplot==1).cl_kutta], 1);
        xFitting = flow.alphaplotval;
        yFitted = polyval(coeffs, xFitting);      
        solutionplot_fit=plot(xFitting, yFitted, 'r', 'LineWidth', 2);
        xlabel('alpha [º]');
        ylabel('cl Kutta-Joukowski');    
        legend('Calculated values', ['Least square fit: Slope ' num2str(coeffs(1)) ', Intercept: ' num2str(coeffs(2))],'Location','southeast');
        grid on
    elseif get(handles.chkcmplot,'Value')==1
    end
end


% --- Executes on button press in chkstreamlines.
function chkstreamlines_Callback(hObject, eventdata, handles)
% hObject    handle to chkstreamlines (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of chkstreamlines


% --- Executes on button press in chkvelocitycontour.
function chkvelocitycontour_Callback(hObject, eventdata, handles)
% hObject    handle to chkvelocitycontour (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of chkvelocitycontour
if get(hObject,'Value')
    set(handles.chkpressure,'Value',0)
end

% --- Executes on button press in tgladdcontour.
function tgladdcontour_Callback(hObject, eventdata, handles)
% hObject    handle to tgladdcontour (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global flow field solution panels solutionplot solutionplot_fit streamplot streamplot_fill contourplot contourplot_fill contourplot_colorbar 
% Delete previous solutions 
try
        delete(streamplot)
        delete(streamplot_fill)
end
try
        delete(contourplot)
        delete(contourplot_fill)
        delete(contourplot_colorbar)
end
try
    delete(solutionplot)
    delete(solutionplot_fit)
end
lectura=get(handles.popsolutions,'Value');
flow.alphacontour=flow.alpha(lectura);
set(handles.txtcontourplot,'String', ['alpha: ' num2str(flow.alphacontour)])
tamany=[str2double(get(handles.lblsizex,'String')),str2double(get(handles.lblsizey,'String'))]; %filed of 5x thickness of profile in x and y direction
N=str2double(get(handles.lblresolution,'String')); %size of field matrix
if flow.type==1
    field=fields_NL(panels,solution(lectura),tamany,N);
elseif flow.type==2
    field=fields_L(panels,solution(lectura),tamany,N);
end
if get(handles.chkvelocitycontour,'Value')==1 
    lectura=get(handles.popxyabs,'Value');
    if lectura==2
        campo=field.u;
        title('X Velocity contours [m/s]')
    elseif lectura==3
        campo=field.v;
        title('Y Velocity contours [m/s]')
    elseif lectura==1
        campo=field.V;
        title('Absolute velocity contours [m/s]')
    end
    [~,contourplot]=contourf(handles.axplot,field.x,field.y,campo,100,'LineStyle','none');
    contourplot_fill=fill(handles.axplot,panels.vertex(:,1),panels.vertex(:,2),'w');
    colormap(handles.axplot,jet);
    contourplot_colorbar=colorbar(handles.axplot);
    axis('image'); % Select the scale of the axes
    xlabel('x [m]');
    ylabel('y [m]');
    delete(legend)
    grid off
end
if get(handles.chkpressure,'Value')==1
    [~,contourplot]=contourf(field.x,field.y,field.p./101325,100,'LineStyle','none');
    contourplot_fill=fill(panels.vertex(:,1),panels.vertex(:,2),'w');
    colormap(handles.axplot,jet);
    contourplot_colorbar=colorbar;
    axis('image'); % Select the scale of the axes
    title('Static pressure [atm]')
    xlabel('x [m]');
    ylabel('y [m]');
    delete(legend)
    grid off
end
if get(handles.chkstreamlines,'Value')==1
    nstreamlines=str2double(get(handles.lblnstreamlines,'String'));
    startx = field.x(1:floor(N/nstreamlines):end,1);
    starty = field.y(1:floor(N/nstreamlines):end,1);
    streamplot=streamline(handles.axplot,field.x,field.y,field.u,field.v,startx,starty);   
    streamplot_fill=fill(panels.vertex(:,1),panels.vertex(:,2),'w');
    title('Streamlines')
    axis('image'); % Select the scale of the axes
    xlabel('x [m]');
    ylabel('y [m]');
    delete(legend);
    grid off
end
if get(handles.chkstreamlines,'Value')&&get(handles.chkvelocitycontour,'Value')
    if lectura==2
        title('X Velocity contours [m/s] and streamlines')
    elseif lectura==3
        title('Y Velocity contours [m/s] and streamlines')
    elseif lectura==1
        title('Absolute velocity contours [m/s] and streamlines')
    end
end
if get(handles.chkpressure,'Value')&&get(handles.chkstreamlines,'Value')
    title('Pressure contours [atm] and streamlines')
end



function lblresolution_Callback(hObject, eventdata, handles)
% hObject    handle to lblresolution (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of lblresolution as text
%        str2double(get(hObject,'String')) returns contents of lblresolution as a double


% --- Executes during object creation, after setting all properties.
function lblresolution_CreateFcn(hObject, eventdata, handles)
% hObject    handle to lblresolution (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function lblxsize_Callback(hObject, eventdata, handles)
% hObject    handle to lblxsize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of lblxsize as text
%        str2double(get(hObject,'String')) returns contents of lblxsize as a double


% --- Executes during object creation, after setting all properties.
function lblxsize_CreateFcn(hObject, eventdata, handles)
% hObject    handle to lblxsize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function lblysize_Callback(hObject, eventdata, handles)
% hObject    handle to lblysize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of lblysize as text
%        str2double(get(hObject,'String')) returns contents of lblysize as a double


% --- Executes during object creation, after setting all properties.
function lblysize_CreateFcn(hObject, eventdata, handles)
% hObject    handle to lblysize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit23_Callback(hObject, eventdata, handles)
% hObject    handle to lblresolution (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of lblresolution as text
%        str2double(get(hObject,'String')) returns contents of lblresolution as a double


% --- Executes during object creation, after setting all properties.
function edit23_CreateFcn(hObject, eventdata, handles)
% hObject    handle to lblresolution (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function lblsizex_Callback(hObject, eventdata, handles)
% hObject    handle to lblsizex (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of lblsizex as text
%        str2double(get(hObject,'String')) returns contents of lblsizex as a double


% --- Executes during object creation, after setting all properties.
function lblsizex_CreateFcn(hObject, eventdata, handles)
% hObject    handle to lblsizex (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function lblsizey_Callback(hObject, eventdata, handles)
% hObject    handle to lblsizey (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of lblsizey as text
%        str2double(get(hObject,'String')) returns contents of lblsizey as a double


% --- Executes during object creation, after setting all properties.
function lblsizey_CreateFcn(hObject, eventdata, handles)
% hObject    handle to lblsizey (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popxyabs.
function popxyabs_Callback(hObject, eventdata, handles)
% hObject    handle to popxyabs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popxyabs contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popxyabs


% --- Executes during object creation, after setting all properties.
function popxyabs_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popxyabs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function lblnstreamlines_Callback(hObject, eventdata, handles)
% hObject    handle to lblnstreamlines (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of lblnstreamlines as text
%        str2double(get(hObject,'String')) returns contents of lblnstreamlines as a double


% --- Executes during object creation, after setting all properties.
function lblnstreamlines_CreateFcn(hObject, eventdata, handles)
% hObject    handle to lblnstreamlines (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in chkpressure.
function chkpressure_Callback(hObject, eventdata, handles)
% hObject    handle to chkpressure (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of chkpressure
if get(hObject,'Value')
    set(handles.chkvelocitycontour,'Value',0)
end


% --- Executes on button press in tglsaveplot.
function tglsaveplot_Callback(hObject, eventdata, handles)
% hObject    handle to tglsaveplot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
filter = {'*.png';'*.epsc';'*.eps'};
[file, path] = uiputfile(filter,'Select a folder, a name and the format to save',['Image' datestr(now,'mm_dd_yy-HH,MM,SS')]);

% Copy axes to new figure and save
outputFileName = [path file];
fignew = figure('Visible','off'); % Invisible figure
newAxes = copyobj(handles.axplot,fignew); % Copy the appropriate axes
set(newAxes,'Position',get(groot,'DefaultAxesPosition')); % The original position is copied too, so adjust it.
set(fignew,'CreateFcn','set(gcbf,''Visible'',''on'')'); % Make it visible upon loading
try
    legend(newAxes,handles.axplot.Legend.String)
end
try
    title(newAxes,handles.axplot.Title.String)
end
try
    colorbar(newAxes)
end
saveas(fignew,outputFileName);
delete(fignew);


% --- Executes on button press in tglsavesolution.
function tglsavesolution_Callback(hObject, eventdata, handles)
% hObject    handle to tglsavesolution (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global solution 
filter = {'*.mat'};
[file, path] = uiputfile(filter,'Select a folder, a name and the format to save',['Solution' datestr(now,'mm_dd_yy-HH,MM,SS')]);

% Copy axes to new figure and save
outputFileName = [path file];
save(outputFileName,'solution');
