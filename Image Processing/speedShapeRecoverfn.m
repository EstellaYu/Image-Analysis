function [velocityPixel]= speedShapeRecoverfn(folderPath)
%% API
% 
% thif function 1) takes input of a folder path with Time Stript Images
%               2) calculate velocity, and recover shape of object
%               3) imput space and time calibration 
%               4) save results to the same folder
%
% *************************************************************************
    
%% Select image stacks
    currentPath = pwd;
    % sample_xt = imread('FIG_XTime_OrthogonalView.tif');
    sample_yt = imread('FIG_TimeY_OrthogonalView.tif'); sample_yt = sample_yt(:,:,1:3);
    cd(folderPath)
    
    opts.Interpreter = 'tex';
    prompt = {'Calibration (\mum/pix): ','frame rate (fps): '};
    dlgtitle = 'Input calibrations';
    dims = [1 30];
    definput = {'2.44','30'};
    answer = inputdlg(prompt,dlgtitle,dims,definput,opts);
    
    % raed the first frame and get the size of image
    xt = imread('XTime_OrthogonalView.tif');
    yt = imread('TimeY_OrthogonalView.tif');
    
    if horizontalDialog == 1
        temp = yt;
        yt = imrotate(xt,90);
        xt = imrotate(temp,90); 
    end
    
%% Velocity Measurement

    h = figure(1);
    set(gcf, 'Unit','Normalized','OuterPosition',[0,0,1,0.9])
    subplot(1,2,1)
    imshow(sample_yt)
    title('Select this region for velocity calculation')
    h1 = subplot(1,2,2);
    imshow(yt)
    title('Select a region (ROI) for slope measurement')
    rect = getrect(h1);
    
    velocity = zeros(1,3);
    h2 = subplot(1,2,1);
    imshow(yt(ceil(rect(2)) : floor(rect(2) + rect(4)) ,...
              ceil(rect(1)) : floor(rect(1) + rect(3)) ))
    i = 1;
    titlestr = sprintf('slope measurement %d (%d/3)\n\n\n\n', i, i);
    title(titlestr)
    
    while i < 4
        [xi,yi] = getline(h2);
        if slopeDialog == 1      
            velocity(i) = str2double(answer{1}) * 1e-6 *...
                str2double(answer{2}) * abs(yi(1) - yi(2))/ abs(xi(1) - xi(2));
            if xi(find(yi == max(yi))) == max(xi)
                velocity(i) = - velocity(i);
            end
            i = i + 1;
            switch i
                case 1
                    titlestr =...
                       sprintf('slope measurement %d (%d/3)\n\n\n\n', i, i);
                case 2
                    titlestr =...
                       sprintf('slope measurement %d (%d/3)\n\nv1 = %.3f (um/s)\n\n', i, i, velocity(1)*1e6);
                case 3
                    titlestr =...
                       sprintf('slope measurement %d (%d/3)\n\nv1 = %.3f (um/s)\nv2 = %.3f (um/s)\n', i, i, velocity(1)*1e6, velocity(2)*1e6);
            end
            title(titlestr)      
        end 
    end
    
    close(h)
    
    velocity = mean(velocity); % [m/s]
    
%% Restore original shape
    %[row, col, # color channels] = size(img)
    [height, width, ~] = size(xt);
    velocityPixel = velocity *1e6 / (str2double(answer{1})) /...
        str2double(answer{2});
    
    newHeight = abs(velocityPixel) * height;
    
    stretchedImage = imresize(xt, [newHeight, width]);
    
    figure (2)
    imshow(stretchedImage)
    
%% Output Result Image

  ei = 0;
  while abs(velocity) < 1
      ei = ei + 1;
      velocity = velocity * 10;
  end
  
  thisfile=sprintf('%s/RestoredShape.tif',folderPath); 
  imwrite(stretchedImage,thisfile);

    msg = sprintf('Calibration: %s (um/pix)\nFrame_rate: %s (fps)\nBubble_Velocity: %.4fE-%d (m/s)\n\n',...
        answer{1}, answer{2}, velocity, ei);
    filename = 'velocity.txt';
    fid = fopen(filename,'wt');
    fprintf(fid, msg);
    fclose(fid);
    
    uiwait(msgbox({'Success! Shape restored image saved as ".tiff"';...
       ''; [msg, 'information above stored as ".txt"']} ,'Success' ,'help'))
    
    cd(currentPath)

%% Dialog functions    
    function i = slopeDialog
        slopeans = questdlg('Measure this slope?', ...
                             'Velocity Measurement', ...
                             'Yes','No','Yes');
        % Handle response
        switch slopeans
            case 'Yes'
                i = 1;
            case 'No'
                i = 0;
        end 
    end

    function horizontal = horizontalDialog
        horizontalans = questdlg('Is the experiment done HORIZONTALLY or VERTICALLY?', ...
                             'Horizontal / Vertical', ...
                             'HORIZONTAL','VERTICAL','VERTICAL');
        % Handle response
        switch horizontalans
            case 'HORIZONTAL'
                horizontal = 1;
            case 'VERTICAL'
                horizontal = 0;
        end 
    end

end