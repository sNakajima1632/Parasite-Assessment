% Shido Nakajima
% Numerical(statistical) analysis of trajectory data. Features average
% speed, mean squared displacement (MSD), MSD with respect to initial
% position, and trajectory step dispersion, along with Pearson correlation
% calculation between the four values and ID of NINV/INV.

% additionally conduct machine learning on average speed, msd, and tsd to
% evaluate if there are any significance in current data for classification.

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

%% average/instantaneous speeds, mean sq displacement, trajectory step dispersion initiation
posX = paraData.x_micron_;
posY = paraData.y_micron_;
instSpeed = paraData.y_micron_;

ParaCategory = {};
ID = {};
AvgSpeed = [];
MSDPrev = [];
MSDOrig = [];
% TSD = square root(mean(step length - average step length)^2)
% espimates the step length dispersion
TSD = [];

% calculate all values for plotting (i=1:114)
for i = 1:length(parasiteidIndex)-1
    % index i integrated into parasiteidIndex for cleanliness
    inow = parasiteidIndex(i);
    ilast = parasiteidIndex(i+1)-1;

   % posXY and speedXY only feature coordinates of current parasiteID
    posXY = [posX(inow:ilast),posY(inow:ilast)];
    speedXY = gradient(posXY')'./gradient(paraData.t_sec_(inow:ilast));

    % list of instantaneous velocities for all points
    instSpeed(inow:ilast) = hypot(speedXY(:,1),speedXY(:,2));
    AvgSpeed = cat(1,AvgSpeed,mean(instSpeed(inow:ilast)));

    % TSD calculation
    stepLength = sqrt(sum(diff(posXY).^2,2));
    TSD = cat(1,TSD,sqrt(mean((stepLength(:)-mean(stepLength)).^2)));
    %{
    stepLength = sqrt(sum(diff(posXY).^2,2));
    TSD = stepLength - mean(StepLength);
    TSD = TSD.^2;
    TSD = mean(TSD);
    TSD = sqrt(TSD);
    %}

    % save parasiteID for saving in new table
    ID = cat(1,ID,paraData.PARASITEID(inow));

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
        ParaCategory = cat(1,ParaCategory,'NINV');
    else
        % msd for reference point being x(t-1) and y(t-1) AKA previous point
        % msd = distance/step
        msdp = mean(sum(diff(posXY).^2,2));
        MSDPrev = cat(1, MSDPrev, msdp);
        % msd for reference point being x(1) and y(1) AKA first point
        % msd = displacement/step
        msdo = mean(sum((posXY(2:end,:)-[posXY(1,1),posXY(1,2)]).^2,2));
        MSDOrig= cat(1, MSDOrig, msdo);
        ParaCategory = cat(1,ParaCategory,'INV');
    end
end

% create table
MSDTable = table(ID,ParaCategory,AvgSpeed,MSDPrev,MSDOrig,TSD);

% create correlation matrix between mean speed, msdprev, msdorig, and tsd
rho = partialcorr([grp2idx(MSDTable.ID) MSDTable.AvgSpeed MSDTable.MSDPrev MSDTable.MSDOrig MSDTable.TSD]);
% remove upper half
rho(logical(triu(ones(size(rho)),1))) = NaN;

%% plot data with interactive selection
f = uifigure('Name','Parasite Analysis','Position',[100 100 1700 800]);

% initiate the plot with first item in paraData
ax = axes(f,'Position',[0.05 0.15 0.3 0.7]);
plot(ax,paraData.t_sec_(currParaIndex:nextParaIndex),instSpeed(currParaIndex:nextParaIndex),'o');
xlabel(ax,'time (sec)');
ylabel(ax,'speed (micron/sec)');
title(ax, 'Parasite speed of ' + string(paraData.PARASITEID(1)));
ylim(ax, [0, max(instSpeed)]);
grid(ax,"on");

% dropdown menu for selecting the parasiteID
paraSelect = uidropdown(f, ...
    'Items',unique(paraData.PARASITEID,'stable'), ...
    'Position',[100 700 75 20], ...
    'ValueChangedFcn',@(src,event) updateParaSpeed(src,parasiteidIndex,ax,paraData));

% average speed plot
axNum = axes(f,'Position',[0.4 0.1 0.125 0.8]);
x=categorical(MSDTable.ParaCategory);
y=MSDTable.AvgSpeed;
swarmchart(axNum,x,y,10);
title(axNum, 'Mean Speed');

% MSD plot
axNum = axes(f,'Position',[0.55 0.1 0.125 0.8]);
y=MSDTable.MSDPrev;
swarmchart(axNum,x,y,10);
title(axNum, 'MSD');

% MSD from initial point plot
axNum = axes(f,'Position',[0.7 0.1 0.125 0.8]);
y=MSDTable.MSDOrig;
swarmchart(axNum,x,y,10);
title(axNum, 'MSD, reference=t(1)');

% TSD from initial point plot
axNum = axes(f,'Position',[0.85 0.1 0.125 0.8]);
y=MSDTable.TSD;
swarmchart(axNum,x,y,10);
title(axNum, 'TSD');

% plot heatmap of correlation coefficient between speed, MSD, and TSD
h = heatmap(rho,...
    'MissingDataColor','k', ...
    'XDisplayLabels',["Category", "Mean Speed", "MSD", "MSD ref=t(1)", "TSD"], ...
    'YDisplayLabels',["Category", "Mean Speed", "MSD", "MSD ref=t(1)", "TSD"]);

%% Use Machine Learning (clustering) to evaluate data
% Hierarchical Clustering evaluation to observe optimal number of clusters
eva = evalclusters(AvgSpeed,'linkage','CalinskiHarabasz','KList',1:114);
disp("Mean speed Hierarchical Clustering result");
disp(eva);
eva = evalclusters(MSDPrev,'linkage','CalinskiHarabasz','KList',1:114);
disp("MSD Hierarchical Clustering result");
disp(eva);
eva = evalclusters(MSDOrig,'linkage','CalinskiHarabasz','KList',1:114);
disp("MSD with reference at t=1 Hierarchical Clustering result");
disp(eva);
eva = evalclusters(TSD,'linkage','CalinskiHarabasz','KList',1:114);
disp("TSD Hierarchical Clustering result");
disp(eva);

% k-means clustering evaluation to observe optimal number of clusters
eva = evalclusters(AvgSpeed,'kmeans','CalinskiHarabasz','KList',1:114);
disp("Mean speed k-Means Clustering result");
disp(eva);
eva = evalclusters(MSDPrev,'kmeans','CalinskiHarabasz','KList',1:114);
disp("MSD k-Means Clustering result");
disp(eva);
eva = evalclusters(MSDOrig,'kmeans','CalinskiHarabasz','KList',1:114);
disp("MSD with reference at t=1 k-Means Clustering result");
disp(eva);
eva = evalclusters(TSD,'kmeans','CalinskiHarabasz','KList',1:114);
disp("TSD k-Means Clustering result");
disp(eva);

%% function updateParaSpeed to update the instantaneous speed plot
function updateParaSpeed(src, parasiteidIndex, ax, paraData)
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