% Shido Nakajima
% Started 9/6/2024
% Assessment from Dr. Ganusov of Texas Biomed

clear;clc;close all;

%% import data from excel as table
paraData = readtable("data-SPZ-in-skin-to-analyze.xlsx");
% sort by movie
paraData = sortrows(paraData,"movie");

% plot movie to visualize. NOT NEEDED
%plot(paraData.movie);

%% list index of where 'movie' value changes
movieIndex = ischange(paraData.movie);
movieIndex = find(movieIndex);

%% list index of where 'PARASITEID' value changes
% loop through parasiteID and extract numbers only for use in ischange()
parasiteidIndex = zeros(length(movieIndex),1);

for i = 1:length(paraData.PARASITEID)
    parasiteidIndex(i) = str2double(extract(paraData.PARASITEID(i), digitsPattern(1,2)));
end

parasiteidIndex = ischange(parasiteidIndex);
% 1 added at beginning for index of first parasiteID
parasiteidIndex = cat(1,1,find(parasiteidIndex));
% last index added for last data index (length+1 for consistency)
parasiteidIndex = cat(1,parasiteidIndex,length(paraData.PARASITEID)+1);

% starting index of current parasite on plot
currParaIndex = 1;

% failed attempt at indexing without loop
%parasiteidIndex = char2double(paraData.PARASITEID);
%parasiteidIndex = find(parasiteidIndex);

%% initiate values for average/instantaneous speeds, mean sq displacement change
% list of xy position
posXY = [paraData.x_micron_(1:parasiteidIndex(2)),paraData.y_micron_(1:parasiteidIndex(2))];
% list of xy velocity
speedXY = gradient(posXY')'./gradient(paraData.t_sec_(1:parasiteidIndex(2)));
% list of instantaneous velocities for all points
instSpeed = hypot(speedXY(:,1),speedXY(:,2));
avgSpeed = mean(instSpeed);
meansqDisp = 0.0;

%% plot positions with interactive selection
f = uifigure('Name','Parasite Position','Position',[100 100 1000 800]);
ax = axes(f,'Position',[0.3 0.3 0.65 0.6]);

% initiate the plot with first item in paraData
plot(ax,paraData.x_micron_(1),paraData.y_micron_(1),'o');
xlabel(ax,'x (micron)');
ylabel(ax,'y (micron)');
title(ax, ['Parasite position at t = ', num2str(paraData.t_sec_(1))]);
xlim(ax, [min(paraData.x_micron_) max(paraData.x_micron_)]);
ylim(ax, [min(paraData.y_micron_) max(paraData.y_micron_)]);
grid(ax,"on");

% UI elements on the figure
% slider for selecting time
paraSld = uicontrol(f,'Style','slider', ...
    'Min',1, 'Max',parasiteidIndex(2)-1, ...
    'Value',1, ...
    'Position',[300 100 550 20], ...
    'SliderStep',[1/(parasiteidIndex(2)-1), 1/(parasiteidIndex(2)-1)]);

% slider title
sldTitle = uicontrol(f,'Style','text', ...
    'Position',[300 125 200 20], ...
    'String','Slider to control time', ...
    'FontSize',12);

% ui panel for selecting the parasiteID and displaying information
paraUI = uipanel(f, ...
    'Title','Parasite Information', ...
    'Position',[20 250 200 500]);

% text information on ui panel paraUI
paraInfo = uicontrol(paraUI,'Style','text', ...
    'Position',[10 5 180 450], ...
    'String',['Parasite ID: ' newline ...
        'Recording:   '+ string(paraData.movie(1)) newline ...
        'Nr:          '+ string(paraData.Nr(1)) newline ...
        'Average Speed: '+ string(avgSpeed) + ' micron' newline ...
        'Instantaneous Speed: ' + string(instSpeed(1)) + ' micron/sec'], ...
    'FontSize',11, ...
    'HorizontalAlignment','left');

% dropdown menu for selecting the parasiteID
paraSelect = uidropdown(paraUI, ...
    'Items',unique(paraData.PARASITEID,'stable'), ...
    'Position',[120 435 75 20], ...
    'ValueChangedFcn',@(src,event) updateParaID(src,paraSld,paraInfo,parasiteidIndex,ax,paraData));

% listener to update plot + info with time slider moving
addlistener(paraSld, 'ContinuousValueChange', @(src, event)updatePlot(src, ax, paraData));
addlistener(paraSld, 'ContinuousValueChange', @(src, event)updateInfo(src, paraInfo, paraData));

%% function updateParaID to update slider + info display + plot + avg/inst speed to match selected parasite upon dropdown selection
function updateParaID(src, sld, infoDisplay, parasiteidIndex, ax, paraData)
    % set current index (EX: 2nd dropdown > ValueIndex=2, currParaIndex = parasiteidIndex(2) = 81)
    index = parasiteidIndex(src.ValueIndex);

    % set slider values (EX: 2nd dropdown > min=81)
    sld.Max = parasiteidIndex(src.ValueIndex + 1)-parasiteidIndex(src.ValueIndex);
    sld.Value = 1;
    sld.SliderStep = [1/(sld.Max-sld.Min), 1/(sld.Max-sld.Min)];

    % calculate new speed and update on global variable
    pXY = [paraData.x_micron_(index:parasiteidIndex(src.ValueIndex+1)),paraData.y_micron_(index:parasiteidIndex(src.ValueIndex+1))];
    sXY = gradient(pXY')'./gradient(paraData.t_sec_(index:parasiteidIndex(src.ValueIndex+1)));
    iSpeed = hypot(sXY(:,1),sXY(:,2));
    aSpeed = mean(iSpeed);

    assignin('base','posXY',pXY);
    assignin('base','speedXY',sXY);
    assignin('base','instSpeed',iSpeed);
    assignin('base','avgSpeed',aSpeed);

    % update global variable currParaIndex
    assignin('base','currParaIndex',index);

    % update plot + info to display first data of the selected parasite
    updatePlot(sld,ax,paraData);
    updateInfo(sld, infoDisplay, paraData)
end

%% function updatePlot to update plot + info with slider moving
function updatePlot(slider, ax, paraData)
    % get value of currParaIndex
    index = evalin('base','currParaIndex');

    % get value that the slider is on
    sliderValue = round(slider.Value);

    % update plot
    cla(ax);
    plot(ax, paraData.x_micron_(index+sliderValue-1),paraData.y_micron_(index+sliderValue-1),'o');
    xlabel(ax,'x (micron)');
    ylabel(ax,'y (micron)');
    title(ax, ['Parasite position at t = ', num2str(paraData.t_sec_(index+sliderValue-1))]);
    xlim(ax, [min(paraData.x_micron_) max(paraData.x_micron_)]);
    ylim(ax, [min(paraData.y_micron_) max(paraData.y_micron_)]);
    grid(ax,"on");
end

%% function updateInfo to update displayed info.
function updateInfo(slider, infoPanel, paraData)
    % get value of currParaIndex,instSpeed,avgSpeed
    index = evalin('base','currParaIndex');
    inst = evalin('base','instSpeed');
    avg = evalin('base','avgSpeed');

    % get value that the slider is on
    sliderValue = round(slider.Value);

    % set paraInfo display
    infoPanel.String = ['Parasite ID: ' newline ...
        'Recording:   '+ string(paraData.movie(index)) newline ...
        'Nr:          '+ string(paraData.Nr(index)) newline ...
        'Average Speed: '+ string(avg) + ' micron' newline ...
        'Instantaneous Speed: ' + string(inst(sliderValue)) + ' micron/sec'];
end