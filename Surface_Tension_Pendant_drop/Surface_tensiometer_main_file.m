%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% wdr_surface_tensiometer_ver1.m
%
% By W. D. Ristenpart
% March 31, 2006
%
% This program calculates the surface tension of a pendant droplet based on
% the algorithm of Rotenberg, Boruvka, and Neumann (1983), as discussed by
% Touhami et al. (1996).  The code is loosely based on earlier work by M.
% Akbarian and E. van Nierop.  Advice and suggestions from E. van Nierop
% are gratefully acknowledged.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Added image process to remove the light reflection due the liquid
% interface
% 
% By Hyoungsoo Kim
% April 30, 2014
% From the typical laboratory setups, sometimes, we can not get the clear
% images. Due to the droplet shape, the light is reflected and reflacted at
% the liquid interface. To minimize this bad signal and to detect the clear
% droplet shape, simple image processing step is added.
%
% Minor correction:
% 1) broken letters are fixed. +- symbol does not appear.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The program performs the following major steps.
%
% 1) Load a user-specified image of a pendant droplet.
% 2) Choose the droplet shape, e.g. typical pendant droplet or inverted
%    droplet. It depends on the density difference.
% 3) Ask the image process. When your droplet has a strong reflection
% inside the droplet, it is quite useful. 
% 4) Ask the user to draw a box around the nozzle (of known diameter), to
%    calibrate the pixel size and specify the orientation. Otherwise, if
%    you know the image resolution. You can also skip this step.
% 5) Perform the calibration procedure. 
% 6) Ask the user to select an area around the pendant droplet
% 7) Performs edge detection on the droplet
% 8) Determines the line of symmetry and obtains two x-y coordinate curves
% 9) Smooths the pixelated x-y curves by removing redundant values
% 10) Calculates the total arc length of each curve
% 11) Converts the x-y coords to r-z coords based on the line of symmetry
% 12) Estimates an initial value Ro for the curvature at (r,z) = (0 0)
% 13) Performs a two-parameter least-squares fitting of the solutions to
%     the 3 coupled 1st order ODEs that describe the drop shape for a given
%     curvature and Bond number, minmizing the variance with the detected r,z
% 14) Asks the user to specify the density difference and nozzle diameter
% 15) Calculates the surface tension and saves the data
%
% The least-squares fitting should take at most a few seconds on a Pentium 4
% or faster machine; any longer indicates problems with the image quality
% or a droplet that is either far from equilibrium or asymmetrical due to
% contamination.
%
% The script requires three functions:
% 1) remove_redundancies.m   (smooths the pixelization of the raw x-y data)
% 2) Fit_surface_tension.m   (calculates drop shape for given Ro and Bond #)
% 3) courbure.m          (describes the three 1st order ODEs for drop shape)
%


% Initialization
clear all; clear all; clear all;
close all; close all;
warning('off','Images:initSize:adjustingMag')
dg = [0 .7 0]; % dark green

format LONGENG

% Default Physical Constants
Drho = 0.06; %0.139; %1.0;  % default value for density difference, g/cm3
nozzleDiam = 819.2; % default value for green nozzle diameter, microns
% nozzleDiam = 1587.5;
g   = 9.81; % m/s^2; change only if you've obtained low earth orbit or higher
B0 = 0.5;  % initial guess for Bond number

%%% Default image resolution
meter_px = 1.e-4;

% Analysis Parameters
cal_minsize = 100;
cal_th = 0;
minsize = 100;     % detected edges smaller than these are discarded
th = 0;             % first try default edge threshold

% Load image
[im_file,path_im] = uigetfile('*.*','Image to be analyzed');
I = imread([path_im,im_file]);

I = im2uint8(I);

% Create new folder to save the result
%% Create the post-processed folder
f1       = fullfile([path_im,im_file(1:end-5),'/']);
if (exist(f1) == 0)
   mkdir (f1);
end

%  if isrgb(I)
%      I = rgb2gray(I);
%  end

if size(I,3) ~= 0
 I = uint8(I(:,:,1));
end
I = imadjust(I);
I_origin = I;
 
image_good = 0;

figure(1), imshow(I);                    

dropshpe = questdlg('The droplet is heavier than the surrounding fluid ?','Determine droplet direction','Yes','No(inverted droplet)','Yes');

%% Image process for the droplet
button_0 = questdlg('Do you want to improve the image for the processing?','Check Image','Yes','No','No');
                    if strcmp(button_0, 'Yes')
                        
                       
                        %% Image process 
                         % Create a gradation background intensity distribution
                            background_left = (I(:,1)+I(:,2))/2;          
                            background_left_min = min(background_left);       background_left_max = max(background_left);

                            background_right = (I(:,end)+I(:,end-1))/2;       
                            background_right_min = min(background_right);     background_right_max = max(background_right);

                            if (background_left_max-background_left_min) >= (background_right_max-background_right_min)
                                background = background_left ;
                            elseif (background_left_max-background_left_min) < (background_right_max-background_right_min)
                                background = background_right ;
                            else
                                return;
                            end

                            for n = 1 : size(I,2)-1;
                                background = horzcat(background, I(:,1));
                            end
                            background_noise = mean(mean(background));

                        I_back     = (background-I);
                        I_temp_0   = imadjust(I_back);
                        I_temp_1   = imfill(I_temp_0,'holes');
                        I_temp_2   = imsharpen(I_temp_1);
                        I_binary   = (mean(max(I_temp_2))*0.3 < I_temp_2);
                        I_binary0  = imfill(I_binary,'holes');
                        
                        I_binary1 = imcomplement(I_binary0);

                        I = im2uint8(I_binary1);         
                        I1 = imcomplement(I);
                        se1    = strel('disk',3);
                        I2   = imerode(I1, se1);
                        I3   = imdilate(I2,se1);
                       
                        
                        I = imcomplement(I3);
                        
                        
                        
                
                        hfig = figure(1);
                        set(hfig, 'Position', [200 200 600 600 ] );
                        subplot(2,2,1);
                        imshow(I_origin);
                        title('Original iamge','FontSize',13) ;
                        subplot(2,2,2);
                        imshow(I);
                        title('Improved image','FontSize',13) ;
                        subplot(2,2,[3 4]);
                        imshow(I_origin-I);
                        title('Subtracted area, i.e. bad signal','FontSize',13) ;
                
                    elseif  strcmp(button_0, 'No')
                        I = I_origin;
                            
                    else 
                        return;
                    end

while image_good == 0;
%%% What is the magnification of the Objective?
    button_1 = questdlg('Do you want to calibrate the image? Otherwise, you can put your calibration information.','Calibration','Image calibration','I know the resolution [pixel/mm]','Image calibration');
    if strcmp(button_1, 'Image calibration')
        
                % Crop nozzle to determine calibration length
  
                hfig = figure(1);
                set(hfig, 'Position', [100 100 1000 700 ] );
                subplot(1,2,1);
                imshow('Needle.bmp');
                title('Needle selection example','FontSize',13) ;
                subplot(1,2,2);
                imshow(I_origin);
                title('Crop image around nozzle for calibration like the left example','FontSize',13) ;
                [Icalib, rectcalib] = imcrop;
                Icalib = I(round(rectcalib(2)):round(rectcalib(2))+round(rectcalib(4)), round(rectcalib(1)):round(rectcalib(1))+round(rectcalib(3)));
                Icalib_2 = medfilt2(imadjust(Icalib)); %increases the contrast of the cropped area
               
                         [bwcalib,cal_th] = edge(Icalib_2,'canny');
                         
%                           bwcalib = bwareaopen(bwcalib,cal_minsize,8);
                    
                    stats = regionprops(bwlabel(bwcalib),'Area','Centroid');
                    centers = [stats.Centroid];
                    xc = centers([1,3]); yc = centers([2,4])
                    nozzle_diam_pixels = sqrt(diff(xc)^2+diff(yc)^2);
                    hcalib = figure(2);
                    imshow(Icalib), hold on;
                    [ycal,xcal] = find(bwcalib==1);
                    plot(xcal,ycal,'.b');
                    line([xc(1) xc(2)],[yc(1) yc(2)],'Color','r');
                    title(['Nozzle Diameter = ', num2str(nozzle_diam_pixels), ' pixels.']);  

                    button_2 = questdlg('Do the detected nozzle edges appear accurate?','Check Image','Yes','No(stop)','direct input','Yes');
                    if strcmp(button_2, 'Yes')
                        image_good = 1;
                        manual_calib = 0;
                    elseif strcmp(button_2, 'No(stop)')
% %                         prompt  = {'Edge Threshold Parameter','Minimum Edge Size'};
% %                         dlgtitle   = 'Image Analysis Parameters';
% %                         lines   = 1;
% %                         def     = {num2str(cal_th),num2str(cal_minsize)};
% %                         answer  = inputdlg(prompt,dlgtitle,lines,def);
% %                         cal_th = eval(answer{1});
% %                         cal_minsize = eval(answer{2});
% %                         close(hcalib);
%                     break;
                    elseif strcmp(button_2,'direct input')
                        manual_calib = 1;
                        close(hcalib);
                        image_good = 1;
                    else            
                        return;
                    end
        
        close (figure(1));
        
        
           % Determine orientation of droplet based on proximity of cropped nozzle to image edges
                    imsz = size(I);
                    Dleft = rectcalib(1);
                    Dup = rectcalib(2);
                    Dright = imsz(2) - (rectcalib(1)+rectcalib(3));
                    Ddown = imsz(1) - (rectcalib(2)+rectcalib(4));
                    Dvals = [Dleft Ddown Dright Dup];
                        
                    % %                             if min(Dvals) == Dup
                    % %                                 I = I';
                    % %                                 I_origin = I_origin';
                    % %                             elseif min(Dvals) == Dright
                    % %                                 I = fliplr(I);
                    % %                                 I_origin = fliplr(I_origin);
                    % %                             elseif min(Dvals) == Ddown
                    % %                                 I = fliplr(I');
                    % %                                 I_origin = fliplr(I_origin');
                    % %                             end

 % Determine orientation of droplet based on proximity of cropped nozzle to image edges
        if strcmp(dropshpe, 'Yes')
         I = rot90(I);
         I_origin = rot90(I_origin);
        elseif strcmp(dropshpe, 'No(inverted droplet)')
         I = rot90(I);
         I = fliplr(I);
         I_origin = rot90(I_origin);
         I_origin = fliplr(I_origin);
        else
            return;
        end

    elseif strcmp(button_1, 'I know the resolution [pixel/mm]')
    %%% Enter the resolution information
    prompt = {'What is the actual resolution [pixel/mm]'};
    dlg_title = 'Resolution';
    num_lines = 1;
    def = {'200'};
    answer = inputdlg(prompt,dlg_title,num_lines,def);
    calibration1 = str2num(answer{1});
                       
        % Determine orientation of droplet based on proximity of cropped nozzle to image edges
        if strcmp(dropshpe, 'Yes')
         I = rot90(I);
         I_origin = rot90(I_origin);
        elseif strcmp(dropshpe, 'No(inverted droplet)')
         I = rot90(I);
         I = fliplr(I);
         I_origin = rot90(I_origin);
         I_origin = fliplr(I_origin);
        else
            return;
        end
        
    manual_calib = 1;
    image_good = 1;     
        
    else            
        return;
    end
    
end

if strcmp(button_2, 'No(stop)')
    h=msgbox('The surface tension calculation is stopped. Please prepare another image.')
%                     break;
end
              
hfig = figure(1);
% Crop nozzle around droplet
figure(hfig); %cla;
set(hfig, 'Position', [100 100 1000 700 ] );
subplot(1,2,1);
imshow('Droplet.bmp');
title('Example: Droplet selection') ;

subplot(1,2,2);
imshow(I_origin);
title(['Crop the droplet like the left example']);
[I2, rect] = imcrop;
I2 = I(round(rect(2)):round(rect(2))+round(rect(4)), round(rect(1)):round(rect(1))+round(rect(3)));
close all;

% clean image
I3 = imadjust(I2);
I4 = medfilt2(I3);

% adjust image analysis parameters until satisfactory
image_good = 0;
while image_good == 0

    % Edge Detection
    if th == 0
        [bw,th]= edge(I4,'log');
    else
        bw = edge(I4,'log',th);
    end
    bw2 = bwareaopen(bw, minsize, 8);

    % Adjust coordinates along line of symmetry
    [y,x] = find(bw2 == 1);
    p = polyfit(x,y,1);
    theta   = atan(p(1));
    d       = abs(y-(p(1).*x+p(2)));  %location of the apex
    apexind     = find(d==min(d));
    xApex = x(apexind);
    yApex = y(apexind);

    % Find r and z values
    m0 = p(1);
    b0 = p(2);
    b1 = y + 1/m0*x;
    x0 = (b1 - b0)/(m0+1/m0);
    y0 = (b1 - b0)/(m0+1/m0)*m0 + b0;
    ind1 = find(y > y0);  % need to modify if orientation not correct
    ind2 = find(y < y0);
    x1 = x(ind1); x2 = x(ind2);
    y1 = y(ind1); y2 = y(ind2);

    % smooth r and z values to remove pixelization
    [xs1,ys1] = Remove_redundancies(x1,y1);
    [xs2,ys2] = Remove_redundancies(x2,y2);
    bs1 = ys1 + 1/m0*xs1;
    xs01 = (bs1 - b0)/(m0+1/m0);
    ys01 = (bs1 - b0)/(m0+1/m0)*m0 + b0;
    bs2 = ys2 + 1/m0*xs2;
    xs02 = (bs2 - b0)/(m0+1/m0);
    ys02 = (bs2 - b0)/(m0+1/m0)*m0 + b0;
    z1 = [sqrt( (xApex - xs01).^2 + (yApex - ys01).^2); 0];
    r1 = [sqrt( (xs1 - xs01).^2 + (ys1 - ys01).^2); 0];
    z2 = [sqrt( (xApex - xs02).^2 + (yApex - ys02).^2); 0];
    r2 = [sqrt( (xs2 - xs02).^2 + (ys2 - ys02).^2); 0];


    % Display detected edges and give user chance to abort
    hfig = figure('Color','w'), imshow(I2), hold on;
    plot([1 rect(3)], polyval(p, [1 rect(3)]),'-r');
    plot(xApex,yApex,'sr')
    plot(x1,y1,'-','Color',dg)
    plot(xs1,ys1,'o','Color',dg);
    plot(x2,y2,'-b',xs2,ys2,'ob');
    % title('Press any key if acceptable or ctrl-c to quit');
    % pause;
%     close all;

    
    button = questdlg('Do the detected droplet edges appear accurate?','Check Image');
    if strcmp(button, 'Yes')
        image_good = 1;
        close all
    elseif strcmp(button, 'No')
        close all
        prompt  = {'Edge Threshold Parameter','Minimum Edge Size'};
        dlgtitle   = 'Image Analysis Parameters';
        lines   = 1;
        def     = {num2str(th),num2str(minsize)};
        answer  = inputdlg(prompt,dlgtitle,lines,def);
        th = eval(answer{1});
        minsize = eval(answer{2});
%         close(hfig);
    else
        return;
    end
    
end


% calculate arclengths
arcL1 = sum(sqrt(diff(xs1).^2 + diff(ys1).^2));
arcL2 = sum(sqrt(diff(xs2).^2 + diff(ys2).^2));

% Determine initial estimate for curvature at apex
maxr1 = max(r1);
closeind = find(abs(r1 - 2/3*maxr1) < 2);
closedists = sqrt( (xApex -x1(closeind)).^2 + (yApex - y1(closeind)).^2);
[closestval, closestind] = min(closedists);
nearapex = [find(xs1 > xs1(closeind(closestind))); length(r1)];
a0 = max(z1(nearapex))/max(r1(nearapex))^2;
R0 = 1 / (2*a0);

% Find least squares fit for curvature and Bond number : Curve 1
c0 = [R0, B0];
opts = optimset('MaxFunEvals',1e4,'TolX',1e-10,'TolFun',1e-10,'MaxIter',500);
[coeffs1,resnorm1] = lsqcurvefit(@Fit_surface_tension, c0, z1, r1,[0 0],[],opts,arcL1);
Ro = coeffs1(1);
B = coeffs1(2);

% Find least squares fit for curvature and Bond number : Curve 2
c0 = coeffs1;
opts = optimset('MaxFunEvals',1e4,'TolX',1e-10,'TolFun',1e-10,'MaxIter',500);
[coeffs2,resnorm2] = lsqcurvefit(@Fit_surface_tension, c0, z2, r2,[0 0],[],opts,arcL2);
Ro(2) = coeffs2(1);
B(2) = coeffs2(2);

% Generate figure
fig1 = figure('Color','w','Units','normalized','Position',[.1 .1 .8 .8]);
ax1 = subplot(2,2,1);
ax2 = subplot(2,2,2); axis off;
ax3 = subplot(2,2,3); box on;
ax4 = subplot(2,2,4); box on;

% Display image and detected curves
axes(ax1); hold on;
imshow(I2);
plot([1 rect(3)], polyval(p, [1 rect(3)]),'-r');
plot(xApex,yApex,'sr')
plot(x1,y1,'-','Color',dg)
plot(xs1,ys1,'o','Color',dg);
plot(x2,y2,'-b',xs2,ys2,'ob');
title(im_file,'Interpreter','none');

% Display curve 1
axes(ax3); hold  on;
plot(z1,r1,'o','Color',dg);
plot(z1,Fit_surface_tension(coeffs1,z1,arcL1),'-','Color','k','LineWidth',2);
title(['Curve 1: Bond = ',num2str(B(1)),', Ro = ', num2str(Ro(1))]);
xlabel('z (Pixels)');
ylabel('r (Pixels)');

% Display curve 2
axes(ax4); hold  on;
plot(z2,r2,'ob');
plot(z2,Fit_surface_tension(coeffs2,z2,arcL2),'-k','LineWidth',2);
title(['Curve 2: Bond = ',num2str(B(2)),', Ro = ', num2str(Ro(2))]);
xlabel('z (Pixels)');
ylabel('r (Pixels)');

% Get user inputs
%manual_calib=1;
% Update the calibration information
if strcmp(button_1, 'Image calibration')

    if manual_calib == 1
    prompt  = {'Name','Density Difference (g/cm^3)','Calibration (pixels/mm)','Temperature (C)','Humidity','Save data? (1=yes, 0=no)'};
else
    prompt  = {'Name','Density Difference (g/cm^3)','Nozzle Diameter (microns)','Temperature (C)','Humidity','Save data? (1=yes, 0=no)'};
end

dlgtitle   = 'Pendant Drop Surface Tensiometer';
lines   = 1;
def     = {'Experiment Name',num2str(Drho),num2str(nozzleDiam),'25','50','1'};
answer  = inputdlg(prompt,dlgtitle,lines,def);
exptname = answer{1};
Drho = eval(answer{2});
nozzleDiam = eval(answer{3});
temperature = eval(answer{4});
humidity = eval(answer{5});
savedata = eval(answer{6});

% calculate calibration ratio
if manual_calib == 1
    calib = nozzleDiam;  % should have been entered in pixels/mm
else
    calib = nozzle_diam_pixels / nozzleDiam * 1000;  % convert measured length into units of pixels / mm
end


  elseif strcmp(button_1, 'I know the resolution [pixel/mm]')
    if manual_calib == 1
        prompt  = {'Name','Density Difference (g/cm^3)','Temperature (C)','Humidity','Save data? (1=yes, 0=no)'};
    else
        prompt  = {'Name','Density Difference (g/cm^3)','Temperature (C)','Humidity','Save data? (1=yes, 0=no)'};
    end

    dlgtitle   = 'Pendant Drop Surface Tensiometer';
    lines   = 1;
    def     = {'Experiment Name',num2str(Drho),'25','50','1'};
    answer  = inputdlg(prompt,dlgtitle,lines,def);
    exptname = answer{1};
    Drho = eval(answer{2});
    nozzleDiam = calibration1;
    temperature = eval(answer{3});
    humidity = eval(answer{4});
    savedata = eval(answer{5});

    % calculate calibration ratio
    if manual_calib == 1
        calib = nozzleDiam;  % should have been entered in pixels/mm
    else
       return;
    end

else
    return;
end


% Calculate dimensional values for curvature and surface tension
Radius_Curvature = Ro / calib;
Bond = B;
gamma = Drho*g*Radius_Curvature.^2 ./ Bond;  % units work out to g / s^2, or mN/m

% Calculate means and standard deviations (note: std only for 2 curves; take with grain of salt)
meangamma = mean(gamma);
stdgamma = std(gamma);
meanRo = mean(Radius_Curvature);
stdRo = std(Radius_Curvature);
meanBond = mean(Bond);
stdBond = std(Bond);

% Display calculated values on figure
pm=char(177);
axes(ax2);
text(0.2,1.0,exptname);
text(0.2,0.9,['\Delta \rho = ', num2str(Drho), ' g/cm^3'],'FontSize',15);
text(0.2,0.8,['T = ', num2str(temperature), ' C'],'FontSize',15);
text(0.2,0.7,['RH = ', num2str(humidity), ' %'],'FontSize',15);
text(0.2,0.6,['cal = ', num2str(calib), ' pixel/mm'],'FontSize',15);
text(0.2,0.5,['Ro = ', num2str(meanRo,'%3.3f'),' ', pm,' ' , num2str(stdRo,'%2.3f'),' mm'],'FontSize',15);
text(0.2,0.4,['Bond = ', num2str(meanBond,'%2.3f'),' ', pm,' ' , num2str(stdBond,'%2.3f'),' (dless)'],'FontSize',15);
text(0.2,0.2,['\gamma = ', num2str(meangamma,'%2.2f'),' ', pm,' ' , num2str(stdgamma,'%2.2f' ),' mN/m'],'FontSize',15)
axis off;

% Save the data and figure if requested
if savedata
    clear('I'); % don't save whole image
%     save([f1,im_file(1:end-4)]);
    saveas(gcf,[f1,im_file(1:end-4),'fig']);
    text(0.2,0.0,'Data saved!')
end
