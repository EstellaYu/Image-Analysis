function []= orthogonalViewfn(folderPath)
%% API
% 
% thif function 1) takes input of a folder path with image sequence
%               2) create folder to store result 
%               3) create 2 orthogonal view images
%
% *************************************************************************
%% Create folder to save the result

    currentPath = pwd;
    
    OVfolderPath = [folderPath,'/FIG_Orthogonal_View'];
    if ~exist(OVfolderPath, 'dir')
       mkdir(OVfolderPath)
    end
    
    
%% Select image stacks

    cd(folderPath)
    prompt = {'Filename Starts with:','extention (e.g.: .tiff, .png): '};
    dlgtitle = 'Input frames to stack up';
    dims = [1 30];
    definput = {'frame_','.tif'};
    answer = inputdlg(prompt,dlgtitle,dims,definput);
    
    files = dir([answer{1},'*',answer{2}]);
    numOfFrames = length(files);
    
    % raed the first frame and get the size of image
    img = imread(files(1).name);
    [height, width, ~] = size(img);
    
    % initialize the timestript image
    col_ts = ceil(width / 2);
    row_ts = ceil(height / 2);
    timeStriptImg = zeros(numOfFrames, width);
    velImg = zeros(height, numOfFrames);
    h = waitbar(0, 'Othogonal View Process Starts!','Name','Processing timestript images ...');
    
    for f = 1:numOfFrames
        filename = files(f).name;
        img = im2double(imread(filename));

        timeStriptImg(f, :) = img(row_ts, :, 1);
        velImg(:,f) = img(:, col_ts, 1);
    
        processmsg = sprintf('Orthogonal View -- current progress: %1.0f%% \n',...
            f/ numOfFrames * 100);
        waitbar(f/ numOfFrames, h, processmsg)
    end
    
    close(h)
    delete(h)
    
%% plot and save time stript images    
    figure(1)
    set(gcf, 'Unit','Normalized','OuterPosition',[0,0,1,0.9])
    subplot(1,2,1)
    imshow(timeStriptImg)
    titlestr = sprintf('X-TIME plot: Stacked width in every frame\n\n');
    xlabel('<----------------    X (width)    ---------------->')
    ylabel('<--------------------------------   Time');
    title(titlestr)
    subplot(1,2,2)
    imshow(velImg)
    titlestr = sprintf('TIME-Y plot: Stacked height in every frame\n\n');
    xlabel('Time    -------------------------------->');
    ylabel('<----------------    Y (height)    ---------------->')
    title(titlestr)
    
    cd(OVfolderPath)
    imwrite(timeStriptImg,'XTime_OrthogonalView.tif'); 
    imwrite(velImg,'TimeY_OrthogonalView.tif');
      
    msg = sprintf('Folder Path:\n%s', OVfolderPath);
    uiwait(msgbox({'Success! Preprocessed Timestript images generated, saved as ".tiff"';...
        ''; msg} ,'Success' ,'help'))
    
    close all
    
    cd(currentPath)
end