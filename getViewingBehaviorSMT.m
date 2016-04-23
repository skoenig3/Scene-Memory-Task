function getViewingBehaviorSMT(FIXATIONFILE,SAMPRATE,imageX,imageY)
% Seth Koenig 06/06/2012 modified on 10/30/2013 for SMT

% Function extracts behavior from eye tracking data from cortex files for free
% viewing of natural scence both novel and manipulated/moved. This part of the
% code is essentially get_Alldata.m. The 2nd part of the function extracts
% angles, distance, and time between fixations.

%Inputs:
%   FIXATIONFILE:   Fixations extracted from cortex e.g. MP120606_1-fixation.mat
%   SAMPRATE:       Sampling rate of eye tracking data

%Outputs:
%   A file with the following name [FIXATIONFILE(1:end-13) '-ViewingBehavior'].
%   File contains the necessary variables and data for succesive functions
%   for modeling viewing behavior using a BCRW.
%   Saved variables are described by variable variablenames.

if nargin < 1
    error('Not enough inputs: function requires FIXATIONFILE');
end
if nargin < 2  %in ms/sampled point
    SAMPRATE = 5;
end
if nargin < 4
    imageX = 800;
    imageY = 600;
end

transitionthreshold = 45; %in degrees for persistance or change in direction
load(FIXATIONFILE);% loads fixaqtionstats

%-----Get fixations and saccade parameters----%
numtrials = length(fixationstats);
fixduration = NaN(numtrials ,75);
distbtwnfix = NaN(numtrials ,75);
timebtwfix = NaN(numtrials,75);
sacduration = NaN(numtrials,75);
sacdist = NaN(numtrials,75);
sacamplitude = NaN(numtrials,75);
sacangle = NaN(numtrials,75);
sacangle_2fix = NaN(numtrials,75);
anglebtwfix = NaN(numtrials,75);
for cndlop = 1:numtrials; 
    fixations = fixationstats{cndlop}.fixations;
    if ~isempty(fixations)
        fixationtimes = fixationstats{cndlop}.fixationtimes;
        saccadetimes =  fixationstats{cndlop}.saccadetimes;
        xy = fixationstats{cndlop}.XY;
        if fixations(1,1) > imageX/2-100 && fixations(1,1) < imageX/2+100 &&...
                fixations(2,1) < imageY/2+100 && fixations(2,1) > imageY/2-100
            fixations(:,1) = [];
            fixationtimes(:,1) = [];
        end
        xy =fixationstats{cndlop}.XY;
        sacduration(cndlop,1:length(saccadetimes)) = diff(saccadetimes,1)'+1;
        fixduration(cndlop,1:length(fixationtimes)) = diff(fixationtimes,1)'+1;
        for i = 1:size(fixations,2)-1
            fixx = round(fixations(1,i)); fixy = round(fixations(2,i));
            fixx(fixx < 1) = 1; fixx(fixx > 800) = 800;
            fixy(fixy < 1) = 1; fixy(fixy > 600) = 600;
            timebtwfix(cndlop,i) = (fixationtimes(2,i+1)+fixationtimes(1,i+1))/2 ...
                -(fixationtimes(2,i)+fixationtimes(1,i))/2;
            if i == 1
                x = fixations(1,i)-400;
                y = fixations(2,i)-300;
                distbtwnfix(cndlop,i) = sqrt(x^2+y^2);
                anglebtwfix(cndlop,i) = atan2(y,x);
                x = fixations(1,i+1)-fixations(1,i);
                y = fixations(2,i+1)-fixations(2,i);
                distbtwnfix(cndlop,i+1) = sqrt(x^2+y^2);
                anglebtwfix(cndlop,i+1) = atan2(y,x);
            else
                x = fixations(1,i+1)-fixations(1,i);
                y = fixations(2,i+1)-fixations(2,i);
                distbtwnfix(cndlop,i+1) = sqrt(x^2+y^2);
                anglebtwfix(cndlop,i+1) = atan2(y,x);
            end
        end
        fixx = round(fixations(1,end)); fixy = round(fixations(2,end));
        fixx(fixx < 1) = 1; fixx(fixx > 800) = 800;
        fixy(fixy < 1) = 1; fixy(fixy > 600) = 600;
        for i = 1:size(saccadetimes,2)
            sacx = xy(1,saccadetimes(1,i):saccadetimes(2,i));
            sacy = xy(2,saccadetimes(1,i):saccadetimes(2,i));
            sacdist(cndlop,i) = sum(sqrt(diff(sacx).^2+diff(sacy).^2));
            sacamplitude(cndlop,i) = sqrt((sacx(1)-sacx(end))^2+(sacy(1)-sacy(end))^2);
            sacangle(cndlop,i) = atan2(diff(sacy(1:2)),diff(sacx(1:2)));
            sacangle_2fix(cndlop,i) =  atan2(diff(sacy(end-1:end)),diff(sacx(end-1:end)));
        end
    end
end

%---Angle between fixations---%
n = (-180:180)*pi/180;
[probanglebtwfix] = hist(anglebtwfix(~isnan(anglebtwfix)),360);
probanglebtwfix = [probanglebtwfix(36:-1:1) probanglebtwfix probanglebtwfix(end:-1:end-36)];
probanglebtwfix = filtfilt(1/18*ones(1,18),1,probanglebtwfix);
probanglebtwfix = probanglebtwfix(37:end-37);
probanglebtwfix = probanglebtwfix/sum(probanglebtwfix);
probanglebtwfix = [probanglebtwfix probanglebtwfix(1)];

%----Direction of Saccade leaving a fixation---%
[probsacangle] = hist(sacangle(~isnan(sacangle)),360);
probsacangle = [probsacangle(36:-1:1) probsacangle probsacangle(end:-1:end-36)];
probsacangle = filtfilt(1/18*ones(1,18),1,probsacangle);
probsacangle = probsacangle(37:end-37);
probsacangle = probsacangle/sum(probsacangle);
probsacangle = [probsacangle probsacangle(1)];

%---Parameter Profiles of Fixations and Saccades---%
variables = {'Dist','Vel','Accel','Rotation'};
fixlen = round(nanmedian(nanmedian(fixduration)))+10;
saclen = round(nanmedian(nanmedian(sacduration)))+10;
numfixes = sum(sum(~isnan(fixduration)));
numsacs = sum(sum(~isnan(sacduration)));
allfixations = zeros(numfixes,fixlen,length(variables));
allsaccades = zeros(numsacs,saclen,length(variables));
persistence.fix = zeros(numfixes,fixlen);
persistence.sac = zeros(numfixes,saclen);
distanceprofile.fix = NaN(numfixes,fixlen);
distanceprofile.sac = NaN(numfixes,saclen);
fixcount = 1;
saccount = 1;
for cndlop = 1:numtrials; 
    fixationtimes = fixationstats{cndlop}.fixationtimes;
    saccadetimes =  fixationstats{cndlop}.saccadetimes;
    xy =fixationstats{cndlop}.XY;
    if ~isempty(fixationtimes)
        if fixationtimes(1,1)<saccadetimes(1,1)
            fixcount2 = fixcount+1;
        else
            fixcount2 = fixcount;
        end
        for i = 1:size(saccadetimes,2);
            x = xy(1,saccadetimes(1,i):saccadetimes(2,i));
            y = xy(2,saccadetimes(1,i):saccadetimes(2,i));
            velx = diff(x);
            vely = diff(y);
            vel = sqrt(velx.^2+vely.^2);
            accel = abs(diff(vel));
            angle = 180*atan2(vely,velx)/pi;
            transitions = abs(diff(angle)) > transitionthreshold;
            vel = vel(1:end-1);
            rot = zeros(1,length(x)-2);
            dist = zeros(1,length(x)-2);
            for a = 1:length(x)-2;
                rot(a) = abs(angle(a)-angle(a+1));
                dist(a) = sqrt((x(a)-x(a+2)).^2 + (y(a)-y(a+2)).^2);
            end
            rot(rot > 180) = rot(rot > 180)-180;
            if  saclen == length(dist);
                timewarp = 1:length(dist);
            else
                timewarp = round(linspace(1,length(dist),saclen));
            end
            if i == size(saccadetimes,2)
                if fixationtimes(1,end) >= saccadetimes(2,end)
                    distanceprofile.sac(fixcount2,:) = vel(timewarp);
                    persistence.sac(fixcount2,:) = transitions(timewarp);
                end
            else
                distanceprofile.sac(fixcount2,:) = vel(timewarp);
                persistence.sac(fixcount2,:) = transitions(timewarp);
            end
            allsaccades(saccount,:,1) = dist(timewarp);
            allsaccades(saccount,:,2) = vel(timewarp);
            allsaccades(saccount,:,3) = accel(timewarp);
            allsaccades(saccount,:,4) = rot(timewarp);
            saccount = saccount+1;
            fixcount2 = fixcount2+1;
        end
        for i = 1:size(fixationtimes,2);
            x = xy(1,fixationtimes(1,i):fixationtimes(2,i));
            y = xy(2,fixationtimes(1,i):fixationtimes(2,i));
            velx = diff(x);
            vely = diff(y);
            vel = sqrt(velx.^2+vely.^2);
            accel = abs(diff(vel));
            angle = 180*atan2(vely,velx)/pi;
            transitions = abs(diff(angle)) > transitionthreshold;
            vel = vel(1:end-1);
            rot = zeros(1,length(x)-2);
            dist = zeros(1,length(x)-2);
            for a = 1:length(x)-2;
                rot(a) = abs(angle(a)-angle(a+1));
                dist(a) = sqrt((x(a)-x(a+2)).^2 + (y(a)-y(a+2)).^2);
            end
            rot(rot > 180) = rot(rot > 180)-180;
            if  fixlen == length(dist);
                timewarp = 1:length(dist);
            else
                timewarp = round(linspace(1,length(dist),fixlen));
            end
            if i == 1
                if fixationtimes(1,1) >= saccadetimes(2,1)
                    distanceprofile.fix(fixcount,:) = vel(timewarp);
                    persistence.fix(fixcount,:) = transitions(timewarp);
                end
            else
                distanceprofile.fix(fixcount,:) = vel(timewarp);
                persistence.fix(fixcount,:) = transitions(timewarp);
            end
            allfixations(fixcount,:,1) = dist(timewarp);
            allfixations(fixcount,:,2) = vel(timewarp);
            allfixations(fixcount,:,3) = accel(timewarp);
            allfixations(fixcount,:,4) = rot(timewarp);
            fixcount = fixcount+1;
        end
    end
end
avgfixation= mean(allfixations,1);
avgfixprofile = zeros(size(avgfixation));
for i = 1:size(avgfixation,3);
    avgfixprofile(:,:,i) = filtfilt(1/3*ones(1,3),1,avgfixation(:,:,i));
    avgfixprofile(:,:,i) = avgfixprofile(:,:,i) - min(avgfixprofile(:,:,i));
    avgfixprofile(:,:,i) = avgfixprofile(:,:,i)/max(avgfixprofile(:,:,i));
end
avgsaccade= mean(allsaccades,1);
avgsacprofile = zeros(size(avgsaccade));
for i = 1:size(avgsaccade,3);
    avgsacprofile(:,:,i) = filtfilt(1/3*ones(1,3),1,avgsaccade(:,:,i));
    avgsacprofile(:,:,i) =  avgsacprofile(:,:,i) - min(avgsacprofile(:,:,i));
    avgsacprofile(:,:,i) = avgsacprofile(:,:,i)/max(avgsacprofile(:,:,i));
end

Statsbyfixation.fixatoinspertrial = sum(~isnan(fixduration),2)';
Statsbyfixation.meanfixationduration = SAMPRATE*nanmean(fixduration);
Statsbyfixation.stdfixationduration = SAMPRATE*nanstd(fixduration);
Statsbyfixation.numfix = sum(~isnan(fixduration));
Statsbyfixation.meansacdistance = nanmean(sacdist);
Statsbyfixation.stdsacdistance = nanstd(sacdist);
Statsbyfixation.numsacs = sum(~isnan(sacdist));

variablenames={
    'fixduration: durations of fixations';
    'distbtwnfix: distance between fixations';
    'timebtwfix: time between fixations';
    'sacduration: saccadedurations';
    'sacdist: saccade are length';
    'sacamplitude: classic saccade amplitude';
    'sacangle: direction of saccade leaving a fixation';
    'sacangle_2fix: angle of saccade entering a fixation';
    'anglebtwfix: angles between successive fixations';
    'n: angles for probability distributions of angles';
    'probanglebtwfix: probability distribution of angles between successive fixations';
    'probsacangle: probability distrubution of direction of saccade leaving a fixation';
    'persistence: persistence of eye moving in previous direction (change < 45 degrees)';
    'variables: variables for fixation and saccade profiles';
    'allfixations: all variables by fixation';
    'allsaccades: all variables by saccade';
    'avgfixprofile: average variable profile of all fixations warped to median fixation duration';
    'avgsacprofile: average variable profile of all saccades warped to median saccade duration';
    'distanceprofile: profile of distances by paired saccaed and fixation';
    'Statsbyfixation: stats for num fixation per trial, fixation duration by fixation#, & sac distance by saccade #';
    };

save([FIXATIONFILE(1:56) 'datafiles\' FIXATIONFILE(end-22:end-13) '-ViewingBehavior.mat'],...
    'fixduration','distbtwnfix','timebtwfix','sacduration','sacdist','sacangle','anglebtwfix',...
    'n','probanglebtwfix','probsacangle','persistence','variables',...
    'allfixations','allsaccades','avgfixprofile','avgsacprofile',...
    'Statsbyfixation','distanceprofile','sacangle_2fix','sacamplitude')