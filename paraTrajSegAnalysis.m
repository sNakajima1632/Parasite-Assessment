% Shido Nakajima
% Trajectory segmentation(splitting) to analyze trajectory dissociated with
% NINV/INV due to the calculated values clustered together in one dataset 
% shows insignificance. 
% To make machine learning meaningful, decided to segment trajectory to create
% more data to feed into machine learning. Also apply clustering and then 
% manual classification if needed, expecting a more obvious difference between datasets.

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