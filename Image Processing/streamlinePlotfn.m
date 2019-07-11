function [] = streamlinePlotfn(SLfolderPath)

currentPath = pwd;

[shiftGuess, up] = upDownShiftGuess;
cd(currentPath)
if up == true
    frontimg = imread('FIG_front_up.tif'); frontimg = frontimg(:,:,1);
    endimg = imread('FIG_end_up.tif'); endimg = endimg(:,:,1);
    frontStreamline = imread('FIG_steamlineFront_up.tif');frontStreamline = frontStreamline(:,:,3);
    endConnection = imread('FIG_backConnect_up.tif');
else
    frontimg = imread('FIG_front_down.tif'); frontimg = frontimg(:,:,1);
    endimg = imread('FIG_end_down.tif'); endimg = endimg(:,:,1);
    frontStreamline = imread('FIG_steamlineFront_down.tif');frontStreamline = frontStreamline(:,:,3);
    endConnection = imread('FIG_backConnect_down.tif'); endConnection = endConnection(:,:,3);
end

if ~exist(SLfolderPath, 'dir')
   mkdir(SLfolderPath)
end

cd(SLfolderPath)
cd ..

%% Select Image Stacks & enter # of shift/ frame

files = fileDialog;
frames = length(files);

figure(1)
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);
subplot(2,2,1)
imshow(frontimg)
title('You should see - front')
subplot(2,2,3)
imshow(endimg)
title('You should see - end')
hold on

[startFrame,endFrame, shiftFrame] = frameselectDialog(frames, shiftGuess);
previewFrontEndShift(startFrame,endFrame, shiftFrame, files)
again = againDialog;

while again == 1
    [startFrame,endFrame, shiftFrame] = frameselectDialog2(frames,startFrame,endFrame, shiftFrame);
    previewFrontEndShift(startFrame,endFrame, shiftFrame, files)
    again = againDialog;
end
hold off
close all 

%% Preview streamline at the front, select # of frames needed to skip for smooth connections

figure(2)
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);
subplot(2,2,1)
imshow(frontStreamline)
title('You should see -- front, adjust entry 1')
subplot(2,2,3)
imshow(endConnection)
title('You should see -- end, adjust entry 1')
hold on 

[connectFrame, stackFrame] = ConnectselectDialog;
[height, width] = previewConnect(startFrame, endFrame, shiftFrame, connectFrame, stackFrame, files,up);
again = againDialog;
while again == 1
    [connectFrame, stackFrame] = ConnectselectDialog2(connectFrame, stackFrame);
    [height, width] = previewConnect(startFrame, endFrame, shiftFrame, connectFrame, stackFrame, files,up);
    again = againDialog;
end
figure(2)
hold off 
close all

%% Processing streamline plot
msgh = waitbar(0, {'Working hard on it...';'';...
    'This might take ~ a minute (or two)'},'Name','Processing...');

thisframe = streamline(startFrame,endFrame, shiftFrame, connectFrame, stackFrame, files, height, width,msgh, up);

close(msgh)

%% save the result streamline plot

save = saveDialog;
if save == 1
    cd(SLfolderPath)
    imwrite(thisframe,'streamline_BubbleFrame.tif');
    
    msg = sprintf('Folder Path:\n%s', SLfolderPath);
    uiwait(msgbox({'Success! Streamline Figure saved as ".tif"';''; msg},'Success' ,'help'))
end

close all

cd(currentPath)

%% functions used in this fn script
    function [shiftGuess, up] = upDownShiftGuess
       if exist('FIG_Orthogonal_View/velocity.txt','file')     
           cd('FIG_Orthogonal_View/')
           fileID = fopen('velocity.txt');
           vel_file = textscan(fileID, '%s %f %s\n');
           fclose(fileID);
           calibration = vel_file{2}(1); vel = vel_file{2}(3);
           shiftGuess = ceil(vel*1e6/calibration);
           up = vel>0;
           cd ..
       else
           up = str2double(questdlg('Object moving UP or DOWN?','Velocity Direction','UP','DOWN','UP'));
           up = logical(strcmp(up, 'UP'));
           if up == true
               shiftGuess = 5;
           else
               shiftGuess = -5;
           end
       end
    end

    function [files, frames] = fileDialog
        prompt = {'Filename Starts with:','extention (e.g.: .tiff, .png): '};
        dlgtitle = 'Input frames to stack up';
        dims = [1 30];
        definput = {'frame_','.tif'};
        answer = inputdlg(prompt,dlgtitle,dims,definput);
        ext = answer{2};
        files = dir([answer{1},'*', ext]);
        frames = length(files);
    end

    function [startFrame, endFrame, shiftFrame] = frameselectDialog(frames, shiftGuess)
        promptEnd = sprintf('End frame (max frame %d):', frames);
        promptShift = sprintf('Object moves # px/frame');
        prompt = {'Start frame:', promptEnd, promptShift};
        dlgtitle = 'Input frames for analyze (start: end)';
        dims = [1 30];
        definput = {'1', num2str(frames), num2str(shiftGuess)};
        answer = inputdlg(prompt,dlgtitle,dims,definput);

        startFrame = str2double(answer{1});
        endFrame = str2double(answer{2});
        shiftFrame = str2double(answer{3});
    end

    function [startFrame, endFrame, shiftFrame] = frameselectDialog2(frames,startFrame, endFrame, shiftFrame)
        promptEnd = sprintf('End frame (max frame %d):', frames);
        promptShift = sprintf('Object moves # px/frame');

        prompt = {'Start frame:', promptEnd, promptShift};

        dlgtitle = 'Input frames for analyze (start: end)';
        dims = [1 30];
        definput = {num2str(startFrame), num2str(endFrame), num2str(shiftFrame)};
        answer = inputdlg(prompt,dlgtitle,dims,definput);

        startFrame = str2double(answer{1});
        endFrame = str2double(answer{2});
        shiftFrame = str2double(answer{3});
    end

    function [connectFrame, stackFrame] = ConnectselectDialog
        promptConnect = sprintf('need to shift # frames to connect the bubble');
        promptStack = sprintf('# of frames to overlap for steamline');
        prompt = {promptConnect, promptStack};
        dlgtitle = 'Streamline Figure pannel connection';
        dims = [1 30];
        definput = {'80', '50'};
        answer = inputdlg(prompt,dlgtitle,dims,definput);

        connectFrame = str2double(answer{1});
        stackFrame = str2double(answer{2});
    end

    function [connectFrame, stackFrame] = ConnectselectDialog2(connectFrame, stackFrame)
        promptConnect = sprintf('need # frames to connect the bubble');
        promptStack = sprintf('# of frames to overlap for steamline');
        prompt = {promptConnect, promptStack};

        dlgtitle = 'Streamline Figure pannel connection';
        dims = [1 30];
        definput = {num2str(connectFrame), num2str(stackFrame)};
        answer = inputdlg(prompt,dlgtitle,dims,definput);

        connectFrame = str2double(answer{1});
        stackFrame = str2double(answer{2});

    end

    function [again] = againDialog
        answer = questdlg('Choose frames again?','Frame Selection','Enter again','Proceed to next step','Enter again');
        again = strcmp(answer, 'Enter again');
    end

    function previewFrontEndShift(startFrame, endFrame, shiftFrame, files)
        preimg = imread(files(startFrame).name); preimg2 = imread(files(startFrame+10).name);
        preimgend = imread(files(endFrame).name); preimgend2 = imread(files(endFrame-10).name);
        
        if shiftFrame > 0
            preimg = preimg(abs(shiftFrame*10) + 1: end, :, 1); preimg2 = preimg2(1: end-abs(shiftFrame*10),:, 1);
            preimgend2 = preimgend2(abs(shiftFrame*10) + 1: end, :, 1); preimgend = preimgend(1: end-abs(shiftFrame*10),:, 1);
        else
            preimg = preimg(1: end-abs(shiftFrame*10),:, 1); preimg2 = preimg2(abs(shiftFrame*10) + 1: end, :, 1);
            preimgend2 = preimgend2(1: end - abs(shiftFrame*10),:, 1); preimgend = preimgend(abs(shiftFrame*10) + 1: end, :, 1);
        end
        
        figure(1)
        subplot(2,2,2)
        imshow(imcomplement(imabsdiff(preimg, preimg2) * 2))
        title('current front')
        subplot(2,2,4)
        imshow(imcomplement(imabsdiff(preimgend, preimgend2) *2))
        title('crrent end')
    end

    function [height, width]= previewConnect(startFrame, endFrame, shiftFrame, connectFrame, stackFrame, files, up)
        
        preImgStart = imread(files(startFrame).name);        
        [height, width, ~] = size(preImgStart);
        
        %required: height > shift * stackFrameNumber + h
        % h < height - shift * stackFrameNumber
        h = height - abs(shiftFrame) * stackFrame; % height of each pannel
        
        concatFrames = [startFrame, startFrame + connectFrame, endFrame - connectFrame, endFrame];
        % initialize image for display, each pannel as a layer
        showImgfront = zeros(h * 2 , width);
        showImgend = zeros(h * 2 , width);
        
        for i = [1, 2]
            if up == true 
                 rect = [0, height - h + 1, width, h - 1];
            else rect = [0, 0, width, h];
            end
            
            preImg1 = imread(files(concatFrames(i)).name);
            preImg1 = imcrop(preImg1, rect); preImg1 = preImg1(:,:,1);
            
            if up == true 
                 showImgfront( h * (i - 1) + 1 : h * i, :) = preImg1;
            else showImgfront( end - h * i + 1 : end - h * (i - 1), :) = preImg1;
            end
            
            for j = 1:2: stackFrame - 1
                if up == true 
                     rect = [0, height - h - 1 - j * shiftFrame, width, h - 1];
                else rect = [0, j * abs(shiftFrame), width, h - 1];
                end
                
                preImg12 = imcrop(imread(files(concatFrames(i)+j).name), rect);
                preImg12 = preImg12(: ,: ,1);
                diffimg = im2double(imcomplement(imabsdiff(preImg1, preImg12)));
                preImg1 = preImg12;
                if up == true
                    if j == 1    
                        showImgfront( h * (i - 1) + 1 : h * i, :) = diffimg;
                    else
                        showImgfront( h * (i - 1) + 1 : h * i, :) =...
                            imabsdiff(showImgfront( h * (i - 1) + 1 : h * i, :), diffimg);
                    end
                else
                    if j == 1    
                        showImgfront(end - h * i + 1: end - h * (i - 1), :) = diffimg;
                    else
                        showImgfront(end - h * i + 1 :...
                                                end - h * (i - 1), :)...
                                    = ...
                                     imabsdiff(showImgfront( end - h * i + 1: ...
                                                end - h * (i - 1), :), diffimg);
                    end
                end
            end
            
            preImg2 = imread(files(concatFrames(i+2)).name);
            preImg2 = imcrop(preImg2, rect);  
            if up == true
                showImgend( h * (i - 1) + 1 : h * i, :) = preImg2(: ,: , 1);
            else 
                showImgend( end - h * i + 1: end - h * (i - 1), :) = preImg2(: ,: , 1); 
            end
        end
        figure(2)
        subplot(2,2,2)
        imshow((showImgfront-min(min(showImgfront)))./(max(max(showImgfront))-min(min(showImgfront))))
        subplot(2,2,4)
        imshow(showImgend/225)
        flip = flipDialog;
        while flip == 1
            showImgfront = 1- showImgfront;
            subplot(2,2,2)
            imshow((showImgfront-min(min(showImgfront)))./...
                (max(max(showImgfront))-min(min(showImgfront))))
            flip = flipDialog;
        end
    end

    function [thisframe]= streamline(startFrame, endFrame, shiftFrame, connectFrame, stackFrame, files, height, width, msgh, up)

        h = height - abs(shiftFrame) * stackFrame; % height of each pannel
       
        % calcualte num of connections n(1) = startFrame
        i = 0;
        concatFrames = startFrame: connectFrame: endFrame + 100 - i;
        while concatFrames(end) + stackFrame > length(files)
            i = i + 1;
            concatFrames = startFrame: connectFrame: (endFrame + 100 - i);
        end
        n = length(concatFrames);
        
        % initialize image for display, each pannel as a layer
        showImg = zeros(h * n, width);
        
        for i = 1: n
            if up == true 
                 rect = [0, height - h - 1, width, h - 1];
            else rect = [0, 0, width, h];
            end
            
            preImg1 = imread(files(concatFrames(i)).name);
            preImg1 = imcrop(preImg1, rect); preImg1 = preImg1(:,:,1);
            
            if up == true
                 showImg( h * (i - 1) + 1 : h * i, :) = preImg1;
            else showImg( end - h * i + 1 : end - h * (i - 1), :) = preImg1;
            end
            
            for j = 1:2: stackFrame - 1
                if up == true
                     rect = [0, height - h - 1 - j * shiftFrame, width, h - 1];
                else rect = [0, j * abs(shiftFrame) + 1, width, h - 1];
                end 
                
                preImg12 = imcrop(imread(files(concatFrames(i)+j).name), rect);
                preImg12 = preImg12(: ,: ,1);
                diffimg = im2double(imcomplement(imabsdiff(preImg1, preImg12)));
                preImg1 = preImg12;
                if up == true
                    if j == 1    
                        showImg( h * (i - 1) + 1 : h * i, :) = diffimg;
                    else
                        showImg( h * (i - 1) + 1 : h * i, :) = imabsdiff(showImg( h * (i - 1) + 1 : h * i, :), diffimg);
                    end
                else
                    if j == 1    
                        showImg( end - h * i + 1: end - h * (i - 1), :) = diffimg;
                    else
                        showImg( end - h * i + 1: end - h * (i - 1), :) =...
                            imabsdiff(showImg( end - h * i + 1: end - h * (i - 1), :), diffimg);
                    end
                end
            end
            
            progressMsg = sprintf('Streamline Plot -- current progress: %1.0f%% \n',...
                i/n * 100);
            waitbar(i/n ,msgh,progressMsg)
            
        end
        
        figure(3)
        set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 0.3, 0.96]);
        imshow((showImg-min(min(showImg)))./(max(max(showImg))-min(min(showImg))))
        flip = flipDialog;
        while flip == 1
            showImg = 1 - showImg;
            imshow((showImg-min(min(showImg)))./(max(max(showImg))-min(min(showImg))))
            flip = flipDialog;
        end
        thisframe = (showImg-min(min(showImg)))./(max(max(showImg))-min(min(showImg)));      
    end

    function [flip] = flipDialog
        answer = questdlg('show complementary image?','Complement Image','(Yes) Flip the color','(NO) Onward >>','(NO) Onward >>');
        flip = strcmp(answer, '(Yes) Flip the color');
    end

    function [save] = saveDialog
        answer = questdlg('Save the Result?','Save Result','YASSSSS','NO','YASSSSS');
        save = strcmp(answer, 'YASSSSS');
    end
end 