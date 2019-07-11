function [] = Movie2Framesfn(moviePath)
%% API 
% this function 1. slices movie into individual frames,
%               2. optionally rotate frames with certain angle 
%                  (with preview)
%               3. stores frames to folder with the same name of the movie

% moviePath = Folder/filename.ext
% startFrame : stepFrame : endFrame = frames slicing 
%
%**************************************************************************
%% Read movie

    vidobj = VideoReader(moviePath);
    folderName = moviePath(1: end - 4);

    % create folder directory if not exist
    if ~exist(folderName, 'dir')
       mkdir(folderName)
    end
    
    frames=vidobj.Numberofframes;
    
%% input start, stp, and end frames to slice

    startFrame = 1; endFrame = frames; stepFrame = 1;
    
    again = 1;
    while again == 1
    [startFrame, stepFrame, endFrame, h, again] = chooseFrames(vidobj,...
                                   frames, startFrame, stepFrame, endFrame);
    end
    
    close all 

for f=startFrame:stepFrame:endFrame 
  thisframe = read(vidobj,f);
  
%% Image rotation 

  if f == startFrame
      % rotation = rotationDialog;
      rotation = true; first = true;
      angle = 0; 
      
      while rotation == true
        rotatedframe = imrotate(thisframe, angle, 'bilinear','crop');
        figure(1)
        imshow(rotatedframe)
        titlestr_rot = sprintf('Start Frame Preview, start = %d', startFrame);
        title(titlestr_rot)
        rotation = rotationDialog(first);
        first = false;
        if rotation == true
            angle = angle + rotationAngle;
        end
      end
      
      thisframe = rotatedframe;
      croppedframe = thisframe;
      
      [height, width, ~ ] = size(thisframe);
      xmin = 0; ymin = 0; rect(3) = height; rect(4) = width;  first = true;
      crop = cropDialog(first);
      
      while crop == true
          rect = getrect;
          xmin = xmin + rect(1);
          ymin = ymin + rect(2);
          croppedframe = imcrop(thisframe, [xmin, ymin, rect(3), rect(4)]);
          figure (1)
          imshow(croppedframe)
        crop = cropDialog(first);
        first = false;
      end
      
      thisframe = croppedframe;
      toall = rotationToAllDialog;
      h = waitbar(0, 'Slicing process starts!','Name','Processing');
  end
  
  if toall == true && f ~= startFrame
      thisframe = imrotate(thisframe, angle,'bilinear','crop');
      thisframe = imcrop(thisframe, [xmin, ymin, rect(3), rect(4)]);
  end
  
  %if mod(f-1,100) == 0
      msg = sprintf('Slicing process -- current progress: %1.0f%% \n',...
          length(startFrame:stepFrame:f)/length(startFrame:stepFrame:endFrame)* 100);
      waitbar(length(startFrame:stepFrame:f)/length(startFrame:stepFrame:endFrame), h, msg)
  %end
  
%% Output Image  
  thisfile=sprintf('%s/frame_%05d.tif',folderName,f); 
  imwrite(thisframe,thisfile);
  
end

close(h)
delete(h)

msg = sprintf('Folder Path:\n%s', folderName);
uiwait(msgbox({'Success! Frames saved as ".tiff"';''; msg} ,'Success' ,'help'))

%% Customer Interactive Functions

    function [again] = againDialog
        answer = questdlg('Choose frames again?',... % question
                          'Frame Selection',... % title
                          '(Yes) Again','(No) Proceed to next step >>',...
                          '(Yes) Again'); % defalt button
        again = strcmp(answer, '(Yes) Again');
    end

    function [startFrame, stepFrame, endFrame, h, again] = ...
            chooseFrames(vidobj, frames, startFrame, stepFrame, endFrame)
        h = figure(1);
        set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0, 1, 0.9]);
        subplot(2,2,1)
        imagesc(read(vidobj,startFrame));
        titlestr = sprintf('Start Frame Preview, start = %d', startFrame);
        title(titlestr,'FontSize',20)
        axis off
        
        subplot(2,2,2)
        imagesc(read(vidobj,startFrame+50));
        titlestr = sprintf('Start + 50 frames');
        title(titlestr,'FontSize',10)
        axis off
        
        subplot(2,2,3)
        imagesc(read(vidobj,endFrame-50));
        titlestr = sprintf('End - 50 frames');
        title(titlestr,'FontSize',10)
        axis off
        
        subplot(2,2,4)
        imagesc(read(vidobj,endFrame))
        titlestr = sprintf('End Frame Preview, end = %d', endFrame);
        title(titlestr,'FontSize',20)
        axis off
        
        again = againDialog;
        if again == 1
            promptend = sprintf('End frame (max frame %d):', frames);
            prompt = {'Start frame:','Step frame: ', promptend};
            dlgtitle = 'Choose frames -- Click "OK" to Preview';
            dims = [1 30];
            definput = {num2str(startFrame),'1', num2str(endFrame)};
            answer = inputdlg(prompt,dlgtitle,dims,definput);

            startFrame = str2double(answer{1});
            stepFrame = str2double(answer{2});
            endFrame = str2double(answer{3});
        end
    end

    function rotation = rotationDialog(first)
        if first == true
            question_r = sprintf('Do you want to rotate the image?');
        else
            question_r = sprintf('Continue to rotate the image?');
        end 
        rotanswer = questdlg(question_r, ...
                             'Image Rotation', ...
                             '(Yes) Rotate','(No) Onward >>','(No) Onward >>');
        % Handle response
        rotation = logical(strcmp(rotanswer, '(Yes) Rotate'));
    end 

    function angle = rotationAngle
        angleanswer = inputdlg('rotation angle (cw):',...
            'Input rotation angle in degree(\^o)',[1 30],{'0.0'});
        angle = -str2double(angleanswer{1});
    end

    function toall = rotationToAllDialog
        toAllanswer = questdlg('Apply to all the frames?', ...
                               'Image Rotation and/or Crop', ...
                               'Yes','No','Yes');
        % Handle response
        switch toAllanswer
            case 'Yes'
                toall = true;
            case 'No'
                toall = false;
        end 
    end 

    function crop = cropDialog(first)
        if first == true
            question_c = sprintf('Do you want to CROP the image?');
        else
            question_c = sprintf('Continue to CROP the image?');
        end 
        cropanswer = questdlg(question_c, ...
                             'Image Cropping', ...
                             '(Yes) Crop','(No) Next Step >>',...
                             '(No) Next Step >>');
        % Handle response
        crop = logical(strcmp(cropanswer, '(Yes) Crop'));
    end


end
