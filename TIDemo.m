function varargout = TIDemo(varargin)



gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @tupaq_OpeningFcn, ...
                   'gui_OutputFcn',  @tupaq_OutputFcn, ...
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






% --- Executes just before tupaq is made visible.
function tupaq_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to tupaq (see VARARGIN)

% Choose default command line output for tupaq
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);
set(handles.figure1, 'units', 'normalized', 'position', [0.05 0.15 0.9 0.8])


global totalFilename;
totalFilename='';
global total;
total=0;
global FilterSize;
global Sigma;
global batchmode;
batchmode=0;
global batchpath;
global result;
result=0;
global loaddirectory;
loaddirectory=0;
iptsetpref('ImshowAxesVisible','on');

% --- Outputs from this function are returned to the command line.
function varargout = tupaq_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;





% --- Executes on button press in pushbutton_startsingle.
function pushbutton_startsingle_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_startsingle (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global totalFilename;
global total;
global EPI;
global specific;
global neuronsTumour;
global neuronsNormal;
global cut_thr;
global do_cut;
global totaltotal;
global batchmode;
global batchpath;
global labeled_total;
global procspecific;
global nr;
global nc;
global TUMOR;
global STR;
global Classification;
global savedirectory;
savedirectory=0;


if isequal(size(total),[1 1])
    errordlg('Open image for testing (File/Open)','File Error');
    return;
end


if ~batchmode
    procspecific=~isequal(size(specific),[1 1]);
end


    str=['processing ' totalFilename ' , please wait...'];

h=waitbar(0.8,str);

%%downscalling
totaltotal=total;
[nr, nc,nn]=size(totaltotal);

  %if get(handles.checkbox24,'Value')==1
          %TargetLoM=str2double(get(handles.edit21,'String'));
          % DScale=4/TargetLoM;
           DScale=str2double(get(handles.edit21,'String'));
           totaltotal= imresize(totaltotal,DScale); 
  %end
 %%%stain normalization
  %%
%  imtool(totaltotal)
% TargetImage = imread('Ref.bmp');
%SCDMatrix = EstUsingSCD( totaltotal );
%Deconvolve( totaltotal, SCDMatrix );
% [ totaltotal ] = Norm(totaltotal, TargetImage, 'SCD');
 
 %imtool(totaltotal)
 
 
 
 
 
 
 
 
 
 
 
 
 
 %%
if get(handles.checkbox16,'Value')==1
    figure(1), set(gcf, 'name', 'Original Image'), imshow(totaltotal), title('Original Image');
  %hFigure16=imtool(total);
  %set(hFigure16,'NumberTitle','off','Name','Original Image');
end
small_total=totaltotal;
%small_total=rgb2gray(total);
 NumOfClasses=str2double(get(handles.edit6,'String'));
  SIGMA=str2double(get(handles.edit7,'String'));

  if get(handles.checkbox12,'Value')==0
     flag=1;
  else
      flag=0;
  end
  if flag==1
      t = cputime;
      [Classification, TUMOR, EPI1,STR] = TumorTesting(neuronsNormal,neuronsTumour,small_total,NumOfClasses,SIGMA,1);
      time = cputime-t

      EPI=small_total;
      EPI(EPI1>0)=255;
             if get(handles.checkbox9,'Value')==1
                  figure(2),  set(gcf, 'name', 'Epithelium'),imshow(EPI1),title('Epithelium');
                  %hFigure9=imtool(EPI);
                 %set(hFigure9,'NumberTitle','off','Name','Epithelium');
              end

              if get(handles.checkbox8,'Value')==1
                  figure(3),  set(gcf, 'name', 'Stroma'), imshow(STR), title('Stroma');
                %hFigure8=imtool(STR);
                %set(hFigure8,'NumberTitle','off','Name','Stroma');
              end
              if get(handles.checkbox17,'Value')==1
                    figure(4),  set(gcf, 'name', 'Tumour Epithelium'), imshow(Classification),title('Tumour Epithelium');
              %   hFigure17=imtool(Classification);
              %  set(hFigure17,'NumberTitle','off','Name','Classificatio result');
              end
              waitbar(0.5,h);
              drawnow;


  elseif isempty(neuronsNormal) || isempty(neuronsTumour)
         h = errordlg('SOMs should be trained');
         %waitbar(0.5);
                     % drawnow;
  elseif (flag == 0)
             
         [Classification, TUMOR, EPI,STR] = TumorTesting(neuronsNormal,neuronsTumour,small_total,NumOfClasses,SIGMA,0);
                      if get(handles.checkbox9,'Value')==1
                          figure(5),  set(gcf, 'name', 'Epithelium'), imshow(EPI),title('Epithelium');
                        % hFigure9=imtool(EPI);
                        % set(hFigure9,'NumberTitle','off','Name','Epithelium');
                       
                      end

                      if get(handles.checkbox8,'Value')==1
                          figure(6),  set(gcf, 'name', 'Stroma'),imshow(STR),title('Stroma');
                     %   hFigure8=imtool(STR);
                     %   set(hFigure8,'NumberTitle','off','Name','Stroma');
                      end
                      if get(handles.checkbox17,'Value')==1
                          figure(7),  set(gcf, 'name', 'Tumour Epithelium'),imshow(Classification),title('Tumour Epithelium');
                        % hFigure17=imtool(Classification);
                        %set(hFigure17,'NumberTitle','off','Name','Classificatio result');
                      end
                      waitbar(0.5,h);
                      drawnow;

   elseif (flag == 2)
             
         [Classification, TUMOR, EPI,STR] = TumorTesting(neuronsNormal,neuronsTumour,small_total,NumOfClasses,SIGMA,2);
                      if get(handles.checkbox9,'Value')==1
                          figure(8),  set(gcf, 'name', 'Epithelium'),imshow(EPI),title('Epithelium');
                        % hFigure9=imtool(EPI);
                        % set(hFigure9,'NumberTitle','off','Name','Epithelium');
                      end

                      if get(handles.checkbox8,'Value')==1
                          figure(9),  set(gcf, 'name', 'Stroma'), imshow(STR),title('Stroma');
                        %hFigure8=imtool(STR);
                        %set(hFigure8,'NumberTitle','off','Name','Stroma');
                      end
                      if get(handles.checkbox17,'Value')==1
                          figure(10),  set(gcf, 'name', 'Tumour Epithelium'),imshow(Classification),title('Tumour Epithelium');
                        % hFigure17=imtool(Classification);
                       % set(hFigure17,'NumberTitle','off','Name','Classificatio result');
                      end
                      waitbar(0.5,h);
                      drawnow;
  end
EPI=EPI1;
        %%
        % FileName11=strcat(totalFilename,'RGB.tif');
        % saveDataName = fullfile(FileName11);
        % imwrite(Classification,saveDataName,'compression','none') ;
         %% Lymphocte romoval. Only works with CD3
         
                 step_size=1;
        %[Istogram] = CreateIstogram(step_size);
        load Istogram;
        %testing
        [Bchannel, ON] = DetectLymph(totaltotal,Istogram, step_size, 0.6);
        %imtool(Bchannel);
      
        
        BBchannel=rgb2gray(Bchannel);
         BBchannel=double(imbinarize(BBchannel));
       
         %%fill holes
         BB=~BBchannel;
     %     figure(11), imshow(BB)
      %  imtool(BB)  
       %  BB=imfill(BB,'holes');
        se = strel('disk',2);
        BB = imclose(BB,se);
     %     figure(111), imshow(BB)
       % imtool(BB)

        
           Classification = Classification.*repmat(uint8(~BB),[1,1,3]);
           Classification(BB>0)=50;
           %%
         %   figure(12), imshow(Classification)
           %imtool(Classification)
  %         TUMOR(BB>0)=0;
         %  figure(13), imshow(TUMOR)
           %imtool(TUMOR)
 %end
  % save tumour epithelium region

   % Classification = Classification.*repmat(uint8(Bchannel),[1,1,3]);
    %Classification(BB>0)=50;
         %FileName11=strcat(totalFilename,'RGB1.tif');
        %saveDataName = fullfile(FileName11);
        %imwrite(Classification,saveDataName,'compression','none') ;
         

if get(handles.checkbox17,'Value')==1
     
end


if get(handles.checkbox12,'Value')==1
% hFigure12=imtool(cut_total);
% set(hFigure12,'NumberTitle','off','Name','Object Clusters Division');
end


 
if ~batchmode
    procspecific=0;
end

close(h);

  if get(hObject,'Value')
        set(handles.pushbutton5,'Enable','On');
      set(handles.pushbutton7,'Enable','On');
  else
      set(handles.pushbutton5,'Enable','off');
      set(handles.pushbutton7,'Enable','off');
      
  end

% --------------------------------------------------------------------
function Untitled_2_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Untitled_1_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function File_Open_Callback(hObject, eventdata, handles)
% hObject    handle to File_Open (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)




% --------------------------------------------------------------------
function File_Open_total_Callback(hObject, eventdata, handles)
% hObject    handle to File_Open_total (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global loaddirectory;

if ~isequal(loaddirectory,0)
    [FileName, PathName] = uigetfile(...
        {'*.jpg;*.tiff;*.tif;*.gif;*.png;*.bmp;*.svs',...
        'Image Files (*.jpg,*.tiff,*.tif,*.gif,*.png,*.bmp;*.svs)';...
        '*.*',  'All Files (*.*)'},'Open image for tumour extraction',loaddirectory);
else
    [FileName, PathName] = uigetfile(...
        {'*.jpg;*.tiff;*.tif;*.gif;*.png;*.bmp;*.svs',...
        'Image Files (*.jpg,*.tiff,*.tif,*.gif,*.png,*.bmp;*.svs)';...
        '*.*',  'All Files (*.*)'},'Open image for tumour extraction');
end

if (~isequal(FileName,0))
    % Read and show image
    h=waitbar(0.30,'Opening Figure...');
    i=imread([PathName FileName]);
    axes(handles.image_total)
    %%


        X = cat(3, i(:,:,1), i(:,:,2), i(:,:,3));
    %%
    imshow(X);
    zoom on;
    close(h);
    
    global totalFilename;
    totalFilename=FileName;
    global total;
    total=X;
    % remember the folder where image was loaded from
    loaddirectory=PathName;
    
end

 %if get(hObject,'Value')
        set(handles.pushbutton5,'Enable','off');
      set(handles.pushbutton7,'Enable','off');
 %end

% --------------------------------------------------------------------
function file_exit_Callback(hObject, eventdata, handles)
% hObject    handle to file_exit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% image pixel indexes off from figures
iptsetpref('ImshowAxesVisible','off');

% exit
delete(handles.figure1)



% --------------------------------------------------------------------
function Batch_starttotal_Callback(hObject, eventdata, handles)
% hObject    handle to Batch_starttotal (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


global batchmode;
global totalFilename;
global total;
global specificFilename;
global specific;
global procspecific;

global batchpath;

% enter batch processing mode
batchmode=1;
% save state of the program
temptotalFilename=totalFilename;
temptotal=total;
if procspecific
    tempspecificFilename=specificFilename;
    tempspecific=specific;
end

% select required folder
totaldir = uigetdir('','Select folder containing images for testing');
if procspecific
    specificdir = uigetdir('','Select folder containing specifically stained images');
end
batchpath = uigetdir('','Select folder where to save result files');

% check whether all folders were selected
   
if procspecific&(isequal(totaldir,0)|isequal(batchpath,0)|isequal(specificdir,0))
    errordlg('Select all three folders','folder error');
elseif (~procspecific)&(isequal(totaldir,0)|isequal(batchpath,0))
    errordlg('Select both folders','folder error');
else
    % The Windows operating system caches small "thumbnails" of the images
    % in folders if "thumbnail" view is selected from the "view" menu in 
    % Windows Explorer. This file must be removed before the analysis.
    if (ispc)
        if (exist([totaldir '\Thumbs.db']))
            delete ([totaldir '\Thumbs.db']);
        end
        if procspecific
            if (exist([specificdir '\Thumbs.db']))
                delete ([specificdir '\Thumbs.db']);
            end
        end
    end
    
    totalfiles=dir(totaldir);
    % folder listing for specifically stained objects to test
    % whether it includes same files as totaldir
    if procspecific
        specificfiles=dir(specificdir);
        specificfiles=specificfiles(3:end);
    end
    % remove folders . and ..
    totalfiles=totalfiles(3:end);
    
    % check that totaldir and specificdir really contain exactly the same
    % filenames
    if procspecific
        if ~isequal([totalfiles.name],[specificfiles.name])
            errordlg('Files in both folders must have the same names.'...
                ,'File Error');
            return;
        end
    end
    
    % process each image in folder
    for iter=1:size(totalfiles,1)
        totalFilename=totalfiles(iter).name;
        if procspecific
            specificFilename=totalfiles(iter).name;
        end
        
        % Open the right image files to variables
        % in PC computers, separator of folders is '\'
        try
            if ispc==1
                total=imread([totaldir '\' totalFilename]);
                if procspecific
                    specific=imread([specificdir '\' totalFilename]);
                end
            else
                % in UNIX, separator is '/'
                total=imread([totaldir '/' totalFilename]);
                if procspecific
                    specific=imread([specificdir '/' totalFilename]);
                end
            end
        catch
            errordlg('The selected folders must contain only images','File error');
            return;
        end

        pushbutton_startsingle_Callback(hObject, eventdata, handles);
    end

    message=['Check folder ' batchpath ' for results'];
    msgbox(message,'Done!'); 
end

% exit batch processing mode
batchmode=0;

% restore state of the program
totalFilename=temptotalFilename;
total=temptotal;
if procspecific
    specificFilename=tempspecificFilename;
    specific=tempspecific;
end
batchmode=0;
procspecific=0;



% --- Executes during object creation, after setting all properties.
function edit6_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu4.
function popupmenu4_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu4 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu4


% --- Executes during object creation, after setting all properties.
function popupmenu4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit7_Callback(hObject, eventdata, handles)
% hObject    handle to edit7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit7 as text
%        str2double(get(hObject,'String')) returns contents of edit7 as a double


% --- Executes during object creation, after setting all properties.
function edit7_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in checkbox8.
function checkbox8_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox8


% --- Executes on button press in checkbox9.
function checkbox9_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox9


% --- Executes on button press in checkbox10.
function checkbox10_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox10


% --- Executes on button press in checkbox12.
function checkbox12_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox12 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox12
  if get(hObject,'Value')
        set(handles.pushbutton3,'Enable','On');
      set(handles.pushbutton4,'Enable','On');
  else
      set(handles.pushbutton3,'Enable','off');
      set(handles.pushbutton4,'Enable','off');
      
  end
if ~get(hObject,'Value')
    set(handles.edit12,'Enable','Off');
    
else
    set(handles.edit12,'Enable','On');
end
if ~get(hObject,'Value')
    set(handles.edit13,'Enable','Off');
    
else
    set(handles.edit13,'Enable','On');
end
if ~get(hObject,'Value')
    set(handles.edit14,'Enable','Off');
    
else
    set(handles.edit14,'Enable','On');
end
if ~get(hObject,'Value')
    set(handles.edit15,'Enable','Off');
    
else
    set(handles.edit15,'Enable','On');
end
if ~get(hObject,'Value')
    set(handles.edit16,'Enable','Off');
    
else
    set(handles.edit16,'Enable','On');
end
if ~get(hObject,'Value')
    set(handles.edit17,'Enable','Off');
    
else
    set(handles.edit17,'Enable','On');
end
% --- Executes on button press in checkbox13.
function checkbox13_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox13 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox13



% --- Executes during object creation, after setting all properties.
function popupmenu5_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




% --- Executes on button press in checkbox14.
function checkbox14_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox14 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox14
if ~get(hObject,'Value')
    set(handles.popupmenu5,'Enable','On');
    
else
    set(handles.popupmenu5,'Enable','Off');
end


% --- Executes on button press in checkbox15.
function checkbox15_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox15 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox15
if ~get(hObject,'Value')
    set(handles.edit6,'Enable','On');
    
else
    set(handles.edit6,'Enable','Off');
end


% --- Executes on button press in checkbox16.
function checkbox16_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox16 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox16


% --- Executes on button press in checkbox17.
function checkbox17_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox17 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox17


% --- Executes on button press in checkbox18.
function checkbox18_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox18 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox18
if ~get(hObject,'Value')
    set(handles.edit7,'Enable','Off');
    
else
    set(handles.edit7,'Enable','On');
end


% --- Executes on button press in checkbox19.
function checkbox19_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox19 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox19



function edit9_Callback(hObject, eventdata, handles)
% hObject    handle to edit9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit9 as text
%        str2double(get(hObject,'String')) returns contents of edit9 as a double


% --- Executes during object creation, after setting all properties.
function edit9_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit12_Callback(hObject, eventdata, handles)
% hObject    handle to edit12 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit12 as text
%        str2double(get(hObject,'String')) returns contents of edit12 as a double


% --- Executes during object creation, after setting all properties.
function edit12_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit12 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit13_Callback(hObject, eventdata, handles)
% hObject    handle to edit13 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit13 as text
%        str2double(get(hObject,'String')) returns contents of edit13 as a double


% --- Executes during object creation, after setting all properties.
function edit13_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit13 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit14_Callback(hObject, eventdata, handles)
% hObject    handle to edit14 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit14 as text
%        str2double(get(hObject,'String')) returns contents of edit14 as a double


% --- Executes during object creation, after setting all properties.
function edit14_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit14 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit15_Callback(hObject, eventdata, handles)
% hObject    handle to edit15 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit15 as text
%        str2double(get(hObject,'String')) returns contents of edit15 as a double


% --- Executes during object creation, after setting all properties.
function edit15_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit15 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit16_Callback(hObject, eventdata, handles)
% hObject    handle to edit16 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit16 as text
%        str2double(get(hObject,'String')) returns contents of edit16 as a double


% --- Executes during object creation, after setting all properties.
function edit16_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit16 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit17_Callback(hObject, eventdata, handles)
% hObject    handle to edit17 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit17 as text
%        str2double(get(hObject,'String')) returns contents of edit17 as a double


% --- Executes during object creation, after setting all properties.
function edit17_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit17 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton3.
function pushbutton3_Callback(hObject, eventdata, handles)

global neuronsNormal;
neuronsNormal=[];
% hObject    handle to pushbutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% mapWidth=str2double(get(handles.edit15,'String'));
%   mapHeight=str2double(get(handles.edit16,'String'));
% startLearningRate=str2double(get(handles.edit17,'String'));
% imagefiles = dir('*.tif');      
% nfiles = length(imagefiles);    % Number of files found
% neuronsNormal=[];
% for ii=1:2:nfiles
%  [J]=NetworkTraining(imagefiles(ii+1).name, imagefiles(ii).name,mapWidth,mapHeight,startLearningRate);
% neuronsNormal=[neuronsNormal; J ];
% end
mapWidth=str2double(get(handles.edit15,'String'));
  mapHeight=str2double(get(handles.edit16,'String'));
startLearningRate=str2double(get(handles.edit17,'String'));
%neuronsNormal=[];
totaldir = uigetdir('','Select folder containing images and their masks for training Normal model');
totalfiles=dir(totaldir);
totalfiles=totalfiles(3:end);
Anorm=[];
for iter=1:2:size(totalfiles,1)
    I1=imread([totaldir '\' totalfiles(iter).name]);
    I2=imread([totaldir '\' totalfiles(iter+1).name]);
    %    I1=imread(totalfiles(iter).name);
	%I2=imread(totalfiles(iter+1).name);
    [J]=NetworkTraining(I2, I1,mapWidth,mapHeight,startLearningRate);
    Anorm=[Anorm; J ];
    %%
%         [r c]=size(Anorm1);
%         Ind1=repmat(1,r,1);
    %%
end
% if (flag == 2)
%     neuronsNormal=[Anorm1 Ind1];
% else
[neuronsNormal]=SOMM(Anorm,mapWidth,mapHeight,10000,startLearningRate);
%save neuronsNormal.mat neuronsNormal
%neuronsNormal
%end
% --- Executes on button press in pushbutton4.
function pushbutton4_Callback(hObject, eventdata, handles)
global neuronsTumour;
neuronsTumour=[];

% hObject    handle to pushbutton4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% mapWidth=str2double(get(handles.edit15,'String'));
%   mapHeight=str2double(get(handles.edit16,'String'));
% startLearningRate=str2double(get(handles.edit17,'String'));
% imagefiles = dir('*.tif');      
% nfiles = length(imagefiles);    % Number of files found
% neuronsNormal=[];
% for ii=1:2:nfiles
%  [J]=NetworkTraining(imagefiles(ii+1).name, imagefiles(ii).name,mapWidth,mapHeight,startLearningRate);
% neuronsNormal=[neuronsNormal; J ];
% end
mapWidth=str2double(get(handles.edit12,'String'));
  mapHeight=str2double(get(handles.edit13,'String'));
startLearningRate=str2double(get(handles.edit14,'String'));
%[neuronsTumour]=NetworkTraining('m3t.bmp', 'm3.bmp',mapWidth,mapHeight,startLearningRate);
%neuronsTumour=[];
totaldir = uigetdir('','Select folder containing images and their masks for training Tumour model');
totalfiles=dir(totaldir);
totalfiles=totalfiles(3:end);
Anormm=[];
for iter=1:2:size(totalfiles,1)
        I11=imread([totaldir '\' totalfiles(iter).name]);
    I22=imread([totaldir '\' totalfiles(iter+1).name]);
    %    I1=imread(totalfiles(iter).name);
	%I2=imread(totalfiles(iter+1).name);
    [JJ]=NetworkTraining(I22, I11,mapWidth,mapHeight,startLearningRate);
    Anormm=[Anormm; JJ ];
    %%
%     [r c]=size(Anorm2);
%     Ind2=repmat(2,r,1);
    %%
end
% if(flag == 2)
%     neuronsTumour=[Anorm2 Ind2];
% else
 [neuronsTumour]=SOMM(Anormm,mapWidth,mapHeight,10000,startLearningRate);
%save neuronsTumour.mat neuronsTumour
% neuronsTumour
%end
set(handles.text28,'String',neuronsTumour);
%set(handles.edit19,'string',neuronsTumour);


% --- Executes on button press in pushbutton5.
function pushbutton5_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global total;
global EPI;
global TUMOR;
global nr;
global nc;
global totaltotal;
global STR;
global Classification;
global totalFilename;

total= imresize(total,[ nr nc]);


   Classification= imresize(Classification,[ nr nc]);
     TUMOR= imresize(TUMOR,[ nr nc]);
  TUMOR=imbinarize(double(TUMOR));
  EPI= imresize(EPI,[ nr nc]);
  EPI=imbinarize(double(EPI));
  STR= imresize(STR,[ nr nc]);
  STR=imbinarize(double(STR));

  
 NoTUMOR=nnz(TUMOR);
 
  ScalSize=str2double(get(handles.edit23,'String'));
 TUMORSize= 1/(ScalSize*ScalSize)*NoTUMOR;

      %  ratio=(number_of_objects/(number_of_objectstotal+number_of_objects))*100;
		message = sprintf('the size of tumor EPI in mm^2 is %2.3f ', TUMORSize);
		msgbox(message); 


% --- Executes on button press in checkbox21.
function checkbox21_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox21 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox21


% --- Executes on button press in checkbox22.
function checkbox22_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox22 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox22



function edit20_Callback(hObject, eventdata, handles)
% hObject    handle to edit20 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit20 as text
%        str2double(get(hObject,'String')) returns contents of edit20 as a double


% --- Executes during object creation, after setting all properties.
function edit20_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit20 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkbox23.
function checkbox23_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox23 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox23
if ~get(hObject,'Value')
    set(handles.edit20,'Enable','On');
    
else
    set(handles.edit20,'Enable','Off');
end



function edit21_Callback(hObject, eventdata, handles)
% hObject    handle to edit21 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit21 as text
%        str2double(get(hObject,'String')) returns contents of edit21 as a double


% --- Executes during object creation, after setting all properties.
function edit21_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit21 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end






function edit22_Callback(hObject, eventdata, handles)
% hObject    handle to edit22 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit22 as text
%        str2double(get(hObject,'String')) returns contents of edit22 as a double


% --- Executes during object creation, after setting all properties.
function edit22_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit22 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu7.
function popupmenu7_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu7 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu7


% --- Executes during object creation, after setting all properties.
function popupmenu7_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkbox25.
% function checkbox25_Callback(hObject, eventdata, handles)
% % hObject    handle to checkbox25 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox25


% --- Executes on button press in pushbutton7.
function pushbutton7_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global total;
global EPI;
global TUMOR;
global nr;
global nc;
global totaltotal;
global STR;
global Classification;
global totalFilename;
global TargetImage;
%% stain normalization

%imtool(totaltotal)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%TargetImage = imread('Ref.bmp');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%save TargetImage.mat TargetImage
load TargetImage.mat
%imtool(TargetImage)
%SCDMatrix = EstUsingSCD( totaltotal );
%Deconvolve( totaltotal, SCDMatrix );
%totaltotal = NormSCD(totaltotal, TargetImage);
totaltotal = NormRGBHist(totaltotal, TargetImage);

%[ totaltotal ] = Norm(totaltotal, TargetImage, 'SCD');
 
%imtool(totaltotal)

%%
total= imresize(totaltotal,[ nr nc]);

   Classification= imresize(Classification,[ nr nc]);
     TUMOR= imresize(TUMOR,[ nr nc]);
  TUMOR=imbinarize(double(TUMOR));
  EPI= imresize(EPI,[ nr nc]);
  EPI=imbinarize(double(EPI));
  STR= imresize(STR,[ nr nc]);
  STR=imbinarize(double(STR));
%%
 totaltotal=total;
  Inew = totaltotal.*repmat(uint8(TUMOR),[1,1,3]);
  Inew(Inew==0)=255;
    %Bnew = totaltotal.*repmat(uint8(~TUMOR),[1,1,3]);
    Bnew = totaltotal.*repmat(uint8(STR),[1,1,3]);
  Bnew(Bnew==0)=255;
  
      if get(handles.checkbox21,'Value')==1
           figure(16),  set(gcf, 'name', 'Tumour Epithelium RGB Component'), imshow(Inew),title('TUMOR EPI RGB component');

      end
  
%%%%%%%%%%
%%%%%%%%%%
  t1 = cputime 

  small_total=rgb2gray(Inew);
        % J=medfilt2(small_total,[15 15]);
        %K=abs(small_total-J);
        %small_total=small_total+K;
      bw=total_prepare(small_total,'otsu');
      clean=total_clean(bw,0);
      
  if get(handles.checkbox23,'Value')==1
           [cut_total, Mer]=total_cutter(clean, small_total,9999, 1); 
  else
           hmin=str2double(get(handles.edit20,'String'));
           [cut_total, Mer]=total_cutter(clean, small_total,hmin, 1); 
  end
      
  %%
    [labeled,N] = bwlabel(cut_total,4);
    tempor = regionprops(labeled,'Area');
    Smallest_size=str2double(get(handles.edit22,'String'));
    idx = find([tempor.Area] > Smallest_size);
    bw = ismember(labeled,idx);
    cut_total=bw>0;
  
  %%
      [labeled_total,number_of_objects] = bwlabel(cut_total);

      labeled_total=labeled_total>0;
      if get(handles.checkbox22,'Value')==1
          figure(18),  set(gcf, 'name', 'Nuclei Cluster Division'), imshow(labeled_total),title('Nuclei cluster division in tumour Epithelium');
         %     hFigure119=imtool(labeled_total);
         %    set(hFigure119,'NumberTitle','off','Name','Nuclei cluster division');
      end
              % FileName11=strcat(totalFilename,'TumourEpiNuclei.tif');
         %saveDataName = fullfile(FileName11);
         %imwrite(labeled_total,saveDataName,'compression','none') ;
 %%
   small_total=rgb2gray(Bnew);
         %J=medfilt2(small_total,[15 15]);
        %K=abs(small_total-J);
        %small_total=small_total+K;
      bw=total_prepare(small_total,'otsu');
      clean=total_clean(bw,0);
      
  if get(handles.checkbox23,'Value')==1
           [cut_total, Mer]=total_cutter(clean, small_total,9999, 1); 
  else
           hmin=str2double(get(handles.edit20,'String'));
           [cut_total, Mer]=total_cutter(clean, small_total,hmin, 1); 
  end
      
      [labeled,N] = bwlabel(cut_total,4);
    tempor = regionprops(labeled,'Area');
    Smallest_size=str2double(get(handles.edit22,'String'));
    idx = find([tempor.Area] > Smallest_size);
    bw = ismember(labeled,idx);
    cut_total=bw>0;
  
  
      [labeled_total,number_of_objectstotal] = bwlabel(cut_total);
      if get(handles.checkbox26,'Value')==1
            figure(19),  set(gcf, 'name', 'Nuclei Cluster Division'), imshow(labeled_total),title('Nuclei cluster division in Stroma');
           %   hFigure1199=imtool(labeled_total);
           %  set(hFigure1199,'NumberTitle','off','Name','Nuclei cluster division');
      end
        %FileName11=strcat(totalFilename,'StromaNuclei.tif');
        % saveDataName = fullfile(FileName11);
        % imwrite(labeled_total,saveDataName,'compression','none') ;
%%
%  the proportion of tissue occupied by tumor

% automated pattern recognition morphometric image analysis to
%quantify histologic tumor and nontumor tissue areas in biospecimen tissue sections

%The percentage of tissue area occupied by tumor varied among patients and tumor types and
%was distributed around medians of
        ratio=(number_of_objects/(number_of_objectstotal+number_of_objects))*100;
		message = sprintf('the proportion of tissue occupied by tumor is %2.2f %%', ratio);
		msgbox(message);
		%return;
        time = cputime-t1





function edit23_Callback(hObject, eventdata, handles)
% hObject    handle to edit23 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit23 as text
%        str2double(get(hObject,'String')) returns contents of edit23 as a double


% --- Executes during object creation, after setting all properties.
function edit23_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit23 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkbox26.
function checkbox26_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox26 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox26


% --- Executes on button press in pushbutton8.
function pushbutton8_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global loaddirectory1;

if ~isequal(loaddirectory1,0)
    [FileName1, PathName1] = uigetfile(...
        {'*.jpg;*.tiff;*.tif;*.gif;*.png;*.bmp;*.svs',...
        'Image Files (*.jpg,*.tiff,*.tif,*.gif,*.png,*.bmp;*.svs)';...
        '*.*',  'All Files (*.*)'},'Open image for tumour extraction',loaddirectory1);
else
    [FileName1, PathName1] = uigetfile(...
        {'*.jpg;*.tiff;*.tif;*.gif;*.png;*.bmp;*.svs',...
        'Image Files (*.jpg,*.tiff,*.tif,*.gif,*.png,*.bmp;*.svs)';...
        '*.*',  'All Files (*.*)'},'Open image for tumour extraction');
end

if (~isequal(FileName1,0))
    % Read and show image
    h=waitbar(0.30,'Opening image...');
    i=imread([PathName1 FileName1]);
    axes(handles.image_total)
    %%


        TargetImage = cat(3, i(:,:,1), i(:,:,2), i(:,:,3));
    %%
    figure(20),  set(gcf, 'name', 'Reference'),imshow(TargetImage),title('Ref for stain Normalization');
    %imshow(TargetImage);
    zoom on;
    close(h);
    %save TargetImage.mat TargetImage

    
end



%global TargetImage;
%% stain normalization

%imtool(totaltotal)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%TargetImage = imread('Ref.bmp');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%save TargetImage.mat TargetImage

% hObject    handle to pushbutton4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% mapWidth=str2double(get(handles.edit15,'String'));
%   mapHeight=str2double(get(handles.edit16,'String'));
% startLearningRate=str2double(get(handles.edit17,'String'));
% imagefiles = dir('*.tif');      
% nfiles = length(imagefiles);    % Number of files found
% neuronsNormal=[];
% for ii=1:2:nfiles
%  [J]=NetworkTraining(imagefiles(ii+1).name, imagefiles(ii).name,mapWidth,mapHeight,startLearningRate);
% neuronsNormal=[neuronsNormal; J ];
% end
%mapWidth=str2double(get(handles.edit12,'String'));
%  mapHeight=str2double(get(handles.edit13,'String'));
%startLearningRate=str2double(get(handles.edit14,'String'));
%[neuronsTumour]=NetworkTraining('m3t.bmp', 'm3.bmp',mapWidth,mapHeight,startLearningRate);
%neuronsTumour=[];
%totaldir = uigetdir('','Select folder containing a single Reference image for stain Normalization');
%totalfiles=dir(totaldir);
%totalfiles=totalfiles(3:end);
%Anormm=[];
%for iter=1:2:size(totalfiles,1)
%        TargetImage=imread([totaldir '\' totalfiles(1).name]);
 %       save TargetImage.mat TargetImage

%    I22=imread([totaldir '\' totalfiles(iter+1).name]);
    %    I1=imread(totalfiles(iter).name);
	%I2=imread(totalfiles(iter+1).name);
%    [JJ]=NetworkTraining(I22, I11,mapWidth,mapHeight,startLearningRate);
%    Anormm=[Anormm; JJ ];
    %%
%     [r c]=size(Anorm2);
%     Ind2=repmat(2,r,1);
    %%
%end
% if(flag == 2)
%     neuronsTumour=[Anorm2 Ind2];
% else
% [neuronsTumour]=SOMM(Anormm,mapWidth,mapHeight,10000,startLearningRate);
%save neuronsTumour.mat neuronsTumour
% neuronsTumour

% --- Executes on button press in checkbox27.
function checkbox27_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox27 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox27
  if get(hObject,'Value')
   
      set(handles.pushbutton8,'Enable','On');
  else
      
      set(handles.pushbutton8,'Enable','off');
      
  end
