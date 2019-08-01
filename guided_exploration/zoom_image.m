% Demo of how to zoom and pan/scroll the image using a imscrollpanel control.
clc;    % Clear the command window.
close all;  % Close all figures (except those of imtool.)
imtool close all;  % Close all imtool figures if you have the Image Processing Toolbox.
clear;  % Erase all existing variables. Or clearvars if you want.
workspace;  % Make sure the workspace panel is showing.
format long g;
format compact;
fontSize = 13;
% IMSHOW has historically had this behavior all the way back to 1993, that is, it shows the whole image. 
% If it can show the whole image at �truesize� or more recently called 100, it does,
% but if not, it shows the whole image.
% We introduced IMTOOL and IMSCROLLPANEL, and now there�s more flexibility with the display
% because the scroll panel shows you that there is more image there beyond what you can see.

% Have you tried using IMSCROLLPANEL in combination with IMSHOW? See example below or in:
% doc imscrollpanel
% 
% It allows you to directly set the magnification, programmatically or through a magnification box.
% When we created IMSCROLLPANEL and the associated magnification controls in IPT5, 
% we hoped that they would satisfy exactly this GUI use case.
% 
% Let me know what you think if you try it, or if you've tried it in the past, why it doesn't meet your needs. 
% We're open to discussing this further with you to see if we can fix the issue or make the better solution easier to discover if it's indeed adequate.
% Jeff Mather - Image Processing Toolbox developer team leader.

% Here's the example from the documentation:
% Create a scroll panel with a Magnification Box and an Overview tool.
hFig = figure('Toolbar', 'none',...
              'Menubar', 'none');
img = imread('/home/akshay/Desktop/CVC/guided_exploration/26.jpg');          
hIm = imshow(img);
          
% hIm = imshow('saturn.png');
hSP = imscrollpanel(hFig,hIm); % Handle to scroll panel.
set(hSP,'Units', 'normalized',...
        'Position', [0, .1, 1, .9])
	
% Add a Magnification Box and an Overview tool.
hMagBox = immagbox(hFig, hIm);
boxPosition = get(hMagBox, 'Position');
set(hMagBox,'Position', [0, 0, boxPosition(3), boxPosition(4)])
imoverview(hIm)

% Get the scroll panel API to programmatically control the view.
api = iptgetapi(hSP);
% Get the current magnification and position.
mag = api.getMagnification();
r = api.getVisibleImageRect();

% Demonstrate scrolling.
% View the top left corner of the image.
message = sprintf('Click OK to view the top left corner of the image');
button = questdlg(message, 'Continue?', 'OK', 'Quit', 'OK');
drawnow;	% Refresh screen to get rid of dialog box remnants.
if strcmpi(button, 'Quit')
   return;
end
api.setVisibleLocation(0.5,0.5)

% Change the magnification to the value that just fits.
message = sprintf('Click OK to change the magnification to the value that just fits');
button = questdlg(message, 'Continue?', 'OK', 'Quit', 'OK');
drawnow;	% Refresh screen to get rid of dialog box remnants.
if strcmpi(button, 'Quit')
   return;
end
api.setMagnification(api.findFitMag())

% Zoom in to 1600% on the dark spot.
message = sprintf('Click OK to zoom in to 1600%% on the dark spot.');
button = questdlg(message, 'Continue?', 'OK', 'Quit', 'OK');
drawnow;	% Refresh screen to get rid of dialog box remnants.
if strcmpi(button, 'Quit')
   return;
end
api.setMagnificationAndCenter(16,306,800)

% Zoom in to 100% with upper left showing.
message = sprintf('Click OK to zoom in to exactly 100%%.');
button = questdlg(message, 'Continue?', 'OK', 'Quit', 'OK');
drawnow;	% Refresh screen to get rid of dialog box remnants.
if strcmpi(button, 'Quit')
   return;
end
api.setMagnificationAndCenter(1,1,100)
