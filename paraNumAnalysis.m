% Shido Nakajima
% Machine Learning Image Analysis using an adjusted version of previously built network for image classification.

clear;clc;close all;

%% import excel data, get index of movie and ID
% same as paraSort.m
paraData = readtable("data-SPZ-in-skin-to-analyze.xlsx");
paraData = sortrows(paraData,"movie");

% list index of where 'PARASITEID' value changes
parasiteidIndex = zeros(length(paraData.PARASITEID),1);
for i = 1:length(paraData.PARASITEID)
    parasiteidIndex(i) = str2double(extract(paraData.PARASITEID(i), digitsPattern(1,2)));
end
parasiteidIndex = ischange(parasiteidIndex);
parasiteidIndex = cat(1,1,find(parasiteidIndex));
parasiteidIndex = cat(1,parasiteidIndex,length(paraData.PARASITEID)+1);
% starting index of current parasite on plot
currParaIndex = 1;
nextParaIndex = parasiteidIndex(2)-1;

%% average/instantaneous speeds, mean sq displacement initiation
posX = paraData.x_micron_;
posY = paraData.y_micron_;
instSpeed = paraData.y_micron_;
AvgSpeed = 1:length(parasiteidIndex)-1;

ID = {};
MSDPrev = [];
MSDOrig = [];

% calculate all values for plotting
for i = 1:length(parasiteidIndex)-1
    % index i integrated into parasiteidIndex for cleanliness
    inow = parasiteidIndex(i);
    ilast = parasiteidIndex(i+1)-1;

   % posXY and speedXY only feature coordinates of current parasiteID
    posXY = [posX(inow:ilast),posY(inow:ilast)];
    speedXY = gradient(posXY')'./gradient(paraData.t_sec_(inow:ilast));

    % list of instantaneous velocities for all points
    instSpeed(inow:ilast) = hypot(speedXY(:,1),speedXY(:,2));
    AvgSpeed(i) = mean(instSpeed(inow:ilast));

    % save MSD by NINV or INV
    id = char(paraData.PARASITEID(inow));
    if (id(1) == 'N')
        % msd for reference point being x(t-1) and y(t-1) AKA previous point
        % msd = distance/step
        msdp = mean(sum(diff(posXY).^2,2));
        MSDPrev = cat(1, MSDPrev, msdp);
        % msd for reference point being x(1) and y(1) AKA first point
        % msd = displacement/step
        msdo = mean(sum((posXY(2:end,:)-[posXY(1,1),posXY(1,2)]).^2,2));
        MSDOrig= cat(1, MSDOrig, msdo);
        ID = cat(1,ID,'NINV');
    else
        % msd for reference point being x(t-1) and y(t-1) AKA previous point
        % msd = distance/step
        msdp = mean(sum(diff(posXY).^2,2));
        MSDPrev = cat(1, MSDPrev, msdp);
        % msd for reference point being x(1) and y(1) AKA first point
        % msd = displacement/step
        msdo = mean(sum((posXY(2:end,:)-[posXY(1,1),posXY(1,2)]).^2,2));
        MSDOrig= cat(1, MSDOrig, msdo);
        ID = cat(1,ID,'INV');
    end
end

% create table
AvgSpeed = AvgSpeed';
MSDTable = table(ID,MSDPrev,MSDOrig,AvgSpeed);

%% plot positions with interactive selection
f = uifigure('Name','Parasite Speed','Position',[100 100 1500 800]);
ax = axes(f,'Position',[0.05 0.15 0.4 0.7]);

% initiate the plot with first item in paraData
plot(ax,paraData.t_sec_(currParaIndex:nextParaIndex),instSpeed(currParaIndex:nextParaIndex),'o');
xlabel(ax,'time (sec)');
ylabel(ax,'speed (micron/sec)');
title(ax, 'Parasite speed of ' + string(paraData.PARASITEID(1)));
ylim(ax, [0, max(instSpeed)]);
grid(ax,"on");

% dropdown menu for selecting the parasiteID
paraSelect = uidropdown(f, ...
    'Items',unique(paraData.PARASITEID,'stable'), ...
    'Position',[100 750 75 20], ...
    'ValueChangedFcn',@(src,event) updateParaID(src,parasiteidIndex,ax,paraData));

% average speed plot
axmsdp = axes(f,'Position',[0.475 0.1 0.15 0.8]);
x=categorical(MSDTable.ID);
y=MSDTable.AvgSpeed;
swarmchart(axmsdp,x,y,10);
title(axmsdp, 'Mean Speed');

% MSD plot
axmsdp = axes(f,'Position',[0.65 0.1 0.15 0.8]);
y=MSDTable.MSDPrev;
swarmchart(axmsdp,x,y,10);
title(axmsdp, 'MSD');

% MSD from initial point plot
axmsdo = axes(f,'Position',[0.825 0.1 0.15 0.8]);
y=MSDTable.MSDOrig;
swarmchart(axmsdo,x,y,10);
title(axmsdo, 'MSD from t(1)');

%% function updateParaID to update slider + info display + plot + avg/inst speed to match selected parasite upon dropdown selection
function updateParaID(src, parasiteidIndex, ax, paraData)
    % set current index (EX: 2nd dropdown > ValueIndex=2, currParaIndex = parasiteidIndex(2) = 81)
    index = parasiteidIndex(src.ValueIndex);
    nextindex = parasiteidIndex(src.ValueIndex+1)-1;
    assignin('base','currParaIndex',index);
    assignin('base','nextParaIndex',nextindex);

    % get values to plot
    iSpeed = evalin('base','instSpeed');

    % update plot
    cla(ax);
    plot(ax, paraData.t_sec_(index:nextindex),iSpeed(index:nextindex),'o');
    xlabel(ax,'time (sec)');
    ylabel(ax,'speed (micron/sec)');
    title(ax, 'Parasite speed of ' + string(paraData.PARASITEID(index)));
    ylim(ax, [0, max(iSpeed)]);
    grid(ax,"on");
end