%% Image Processing API 
% open folder ==> obtain path
% slice .mov to .tiff files
% rotate if needed
%
%*************************************************************************

%% FUNCTION 
close all 
clc
clear all


%% welcoming pop-up

uiwait(msgbox({'Hi! Welcome to Image Processing (v.0)!';...
                '';...
                'I can help to: ';...
                '';...
                '1. slice ".mov" file and save individual frames';...
                '   1) slice';... % select directory, slice '.mov'
                '   2) rotate (optional)';... % perview available
                '   3) crop (optional)';...
                '   4) save each frame as ".tif"';...
                '';...
                '2. orthogonal view';...
                '   1) loop through each frame, produce Time Stript Images';... % select directory
                '   2) read Time Stript Images';... 
                '   3) calibrate  & calcualte velocity';...
                '   4) recover object shape';...
                '   5) save shape as ".tif", velocity as ".txt"';...
                '   ';...
                '3. stremline plot';...
                '   1) select start and end frames';... % select directory
                '   2) input shift frame # to accomodate object velocity';... 
                '   3) select connection frame number, and # of frames to stack for streamline';...
                '   4) save streamline as ".tif"';...
                '   '...
                },'Image Processing','help'));
            
            
%% 1. slicing movie
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1) select directory to movie path
% 2) rotate, optional (preview availalbe)
% 3) crop, optional 
% 4) save to same directory (will automatically create a folder if needed)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sliceAnswer = questdlg('1. Process Movie Slicing?', ...
                     'Slice Movie to Frames',...
                     'Yes','No','Yes');

if strcmp(sliceAnswer,'Yes') == 1 
    [file,path] = uigetfile('*.mov','Selec a movie file');
    moviePath = strcat(path,file);
    Movie2Framesfn(moviePath)
end % end slicing movie

close all

%% 2. Orthogonal View
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1) select directory to image sequence
% 2) find image with "start-with-name" and extension (.tif, .png, ect)
% 3) produce orthogonal view images 
% 4) save to new directory FIG_Orthogonal_View 
%    (will automatically create a folder if needed)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if strcmp(sliceAnswer,'Yes') == 1
    OVanswer = questdlg('2. Continue to Process Orthogonal View?', ...
                     'Orthogonal View',...
                     'Yes','No','No');
else
    OVanswer = questdlg('2. Process Orthogonal View?', ...
                     'Orthogonal View',...
                     'Yes','No','Yes');
end


if strcmp(OVanswer,'Yes') == 1
    
    % prompt if there is preprocessed Time Stript Images
    S.Default = 'Yes'; 
    S.IconString = 'custom'; 
    S.IconData = imread('FIG_sampe_TimeStript.tif'); 
    TSIanswer=buttondlg({'2. 1) Do you have preprocessed Time Stipt Images?';...
                         '';...
                         '      Yes  ==>  select folders containing the TSI';...
                         '      No   ==>  select folders of sliced frames to produce TSI'},...
                         'Question','Yes','No',S);
    
    if strcmp(TSIanswer, 'Yes') == 1
        if exist('moviePath','var')
            OVfolderPath = [moviePath(1: end -4),'/FIG_Orthogonal_View'];
        else
            OVfolderPath = uigetdir('Selec folder for Time Stript Images');
        end
        speedShapeRecoverfn(OVfolderPath)
    else
        if exist('moviePath','var')
            OVfolderPath = moviePath(1: end -4);
        else
            OVfolderPath = uigetdir('Selec folder for image sequence');
        end
        
        orthogonalViewfn(OVfolderPath)
        OVfolderPath = [OVfolderPath, '/FIG_Orthogonal_View'];
        speedShapeRecoverfn(OVfolderPath)
    end
end % end orthogonal view processing 

close all

%% 3. Streamline Plot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1) select directory to image sequence
% 2) find image with "start-with-name" and extension (.tif, .png, ect)
% 3) produce streamline plot
% 4) save to new directory FIG_Streamline_Plot 
%    (will automatically create a folder if needed)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
if (strcmp(sliceAnswer,'Yes') == 1 || strcmp(OVanswer,'Yes') == 1)
    SLanswer = questdlg('3. Continue to Process Streamline Plot?', ...
                     'Streamline Plot',...
                     'Yes','No','No');
else
    SLanswer = questdlg('3. Process Streamline Plot?', ...
                     'Streamline Plot',...
                     'Yes','No','Yes');
end

if strcmp(SLanswer, 'Yes') == 1
    if exist('moviePath','var')
        SLfolderPath = [moviePath(1: end - 4),'/FIG_Streamline_Plot'];
    elseif exist('OVfolderPath','var')
        SLfolderPath = [OVfolderPath(1: end - 15),'Streamline_Plot'];
    else
        SLfolderPath = uigetdir('Selec folder for Time Stript Images');
        SLfolderPath = [SLfolderPath,'/FIG_Streamline_Plot'];
    end
    
    streamlinePlotfn(SLfolderPath)
end

