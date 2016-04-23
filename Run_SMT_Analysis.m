%% Runs whole Scene Memory Task (SMT) analysis
% Data is arranged in a VERY SPECIFIC manner. Replace my directories with
% your directories by searching "C:\Users\seth.koenig\Documents\". First
% Create the following folder "Matlab\Scene Memory Task" folder to place following
% Matlab functions and folders...
%
% The following are functions necessary to run all sections of code
%
% 1. getSMTeyedat.m: grab, calibrate, and analysze scan paths from cortex files
% 2. Cluster_Fixation.m: detects fixations and saccades
% 3. getSalienceMap: algorithm creates salience maps
% 4. fixationSalience_and_significanceSMT.m: calcualtes salience at fixations
% 5. getViewingBehaviorSMT: calculates viewing behavior statstics from fixation files
%
% The following folders should be created before analysis
%
% 1. "ITMCNDfiles": Place respecitive cortex .ITM and .CND files here
% 2. "Image Sets": Place image set folders containing images
% 3. "datafiles": leave empty for now
% 4. "cortexfiles": leave empty for now
% 5. "Figures": will store figures here
%
% The following are major sections of the codes. See comments in sections for more details.
% Searching [#] can be used to get to that section. Sections can be run
% independently if but have to have been run in order first.
%
% [0] Import Data and Image files
% [1] Get Fixations and Saccades from Cortex Files
% [2] Get the salience maps for multiple task sets
% [2.5] Get Average Saliency Map
% [3] Calculate Salience at Each Fixation
% [4] Calculate Average Salience at Each Fixation Across Multiple Sets and Monkeys
% [5] Calculate Similarity in Fixation Locations with KL Divergence
% [6] Calculate Viewing Behavior Statistics (e.g. fixation duration and saccade amplitude)
% [7] Combine Viewing Behavior Across Monkeys and Lesion Status

%% Should be run every time if new data to grab files from network and get monkey names and lesion status
%---[0] Import Data and Image files---%
if exist('R:\Buffalo Lab\eblab\Cortex Programs\Scene Memory Task\'); %check if network is connected
    % grab which files are for which set, monkey, and lesion status
    sourcefile ='R:\Buffalo Lab\eblab\Cortex Programs\Scene Memory Task\Sets that Monkeys ran on.xls';
    destinationfile = 'C:\Users\seth.koenig\Documents\MATLAB\Scene Memory Task\Sets that Monkeys ran on.xls';
    sd = dir(sourcefile); %get date file was modified
    dd = dir(destinationfile); %get date file was modified
    
    if ~all(sd.date == dd.date) %if dates are not the same then copy file
        copyfile(sourcefile,destinationfile)
    end
    
    % copy cortex data files from network if not already on your computer
    [~,~,raw] = xlsread(destinationfile);
    for row = 2:size(raw,1);
        if ~strcmpi(raw(row,2),'i'); %i here means ignore file
            cortexfile = raw{row,1};
            if ~exist(['C:\Users\seth.koenig\Documents\MATLAB\Scene Memory Task\cortexfiles\' cortexfile])
                if strcmpi(cortexfile(1:2),'PW')
                    source = ['R:\Buffalo Lab\Cortex Data\Vivian\' cortexfile];
                    destination = ['C:\Users\seth.koenig\Documents\MATLAB\Scene Memory Task\cortexfiles\' cortexfile];
                    copyfile(source,destination);
                elseif strcmpi(cortexfile(1:2),'TT')
                    source = ['R:\Buffalo Lab\Cortex Data\Timmy\' cortexfile];
                    destination = ['C:\Users\seth.koenig\Documents\MATLAB\Scene Memory Task\cortexfiles\' cortexfile];
                    copyfile(source,destination);
                elseif strcmpi(cortexfile(1:2),'MP')
                    source = ['R:\Buffalo Lab\Cortex Data\Peepers\' cortexfile];
                    destination = ['C:\Users\seth.koenig\Documents\MATLAB\Scene Memory Task\cortexfiles\' cortexfile];
                    copyfile(source,destination);
                elseif strcmpi(cortexfile(1:2),'JN')
                    source = ['R:\Buffalo Lab\Cortex Data\Guiseppe\' cortexfile];
                    destination = ['C:\Users\seth.koenig\Documents\MATLAB\Scene Memory Task\cortexfiles\' cortexfile];
                    copyfile(source,destination);
                elseif strcmpi(cortexfile(1:2),'WR')
                    source = ['R:\Buffalo Lab\Cortex Data\Wilbur\' cortexfile];
                    destination = ['C:\Users\seth.koenig\Documents\MATLAB\Scene Memory Task\cortexfiles\' cortexfile];
                    copyfile(source,destination);
                elseif strcmpi(cortexfile(1:2),'RR')
                    source = ['R:\Buffalo Lab\Cortex Data\Red\' cortexfile];
                    destination = ['C:\Users\seth.koenig\Documents\MATLAB\Scene Memory Task\cortexfiles\' cortexfile];
                    copyfile(source,destination);
                elseif strcmpi(cortexfile(1:2),'BZ')
                    source = ['R:\Buffalo Lab\Cortex Data\Bizzie\' cortexfile];
                    destination = ['C:\Users\seth.koenig\Documents\MATLAB\Scene Memory Task\cortexfiles\' cortexfile];
                    copyfile(source,destination);
                else
                    disp('could not find cortexfile. Check source folder and file name')
                end
            end
        end
    end
end

%now create easily accessible data structure to understand monkeys, image sets, and lesion status
%code runs independent of code above if not network connection.
[~,~,raw] = xlsread('C:\Users\seth.koenig\Documents\MATLAB\Scene Memory Task\Sets that Monkeys ran on.xls');
allcortexfiles = {};
rowcount = 1;
for row = 2:size(raw,1);
    if ~strcmpi(raw(row,2),'i'); %i here means ignore file
        allcortexfiles{rowcount,1} = raw{row,1};
        allcortexfiles{rowcount,2} = raw{row,2};
        allcortexfiles{rowcount,3} = raw{row,3};
        rowcount = rowcount+1;
    end
end
save('C:\Users\seth.koenig\Documents\MATLAB\Scene Memory Task\allcortexfiles.mat','allcortexfiles')
clar
%%
%---[1] Get Fixations and Saccades from Cortex Files---%

SMT_dir = 'C:\Users\seth.koenig\Documents\MATLAB\Scene Memory Task\'; %parent directory
cd(SMT_dir)

load('allcortexfiles.mat'); %easily accessible data structure from above

imageX = 800; %horizontal pixel size of images
imageY = 600; %vertical pixel size of images

for crow = 1:size(allcortexfiles,1);
    cortexfile = allcortexfiles{crow,1};
    imset = allcortexfiles{crow,3};
    if imset < 10
        image_dir = [SMT_dir  'Image Sets\SMT00' num2str(imset)];
        CNDFile = [SMT_dir 'ITMCNDfiles\SMT00' num2str(imset) '.CND'];
        ITMFile = [SMT_dir 'ITMCNDfiles\SMT00' num2str(imset) '.ITM'];
    else
        image_dir = [SMT_dir  'Image Sets\SMT0' num2str(imset)];
        CNDFile = [SMT_dir 'ITMCNDfiles\SMT0' num2str(imset) '.CND'];
        ITMFile = [SMT_dir 'ITMCNDfiles\SMT0' num2str(imset) '.ITM'];
    end
    if strcmpi(cortexfile,'PW131024.1')
        CNDFile = [SMT_dir 'ITMCNDfiles\PWSMT00' num2str(imset) '.CND']; %error compiling cortex smt.sav file so had to rewrite for this trial
    end
    getSMTeyedat([SMT_dir 'cortexfiles\' cortexfile],ITMFile,CNDFile,imageX,imageY)
    % grab eye data, calibrate, and detect fixations and saccades
end

%%
%---[2] Get the salience maps for multiple task sets---%
SMT_dir = 'C:\Users\seth.koenig\Documents\MATLAB\Scene Memory Task\'; %parent directory
for imset = 1;
    if imset < 10
        image_dir = [SMT_dir  'Image Sets\SMT00' num2str(imset)];
    else
        image_dir = [SMT_dir  'Image Sets\SMT0' num2str(imset)];
    end
    cd(image_dir)
    dirData = dir(image_dir);
    dirIndex = [dirData.isdir];
    fileList = {dirData(~dirIndex).name}';
    imageindexes = [];
    for i = 1:length(fileList)
        bmps = strfind(fileList{i},'bmp');
        if ~isempty(bmps)
            if double(fileList{i}(bmps-2)) <= 57 %ascii for number
                imageindexes = [imageindexes i];
            end
        end
    end
    for i = 1:length(imageindexes)
        imagefile = fileList{imageindexes(i)};
        getSalienceMap(imagefile) %creates and saves salience maps in Image Sets folders
    end
end
%%
%---[2.5] Get Average Saliency Map---%
SMT_dir = 'C:\Users\seth.koenig\Documents\MATLAB\Scene Memory Task\'; %parent directory

imageX = 800; %horizontal pixel size of images
imageY = 600; %vertical pixel size of images
avgmap = zeros(imageY,imageX);
for imset = 1:5;
    if imset < 10
        image_dir = [SMT_dir  'Image Sets\SMT00' num2str(imset)];
    else
        image_dir = [SMT_dir  'Image Sets\SMT0' num2str(imset)];
    end
    cd(image_dir)
    dirData = dir(image_dir);
    dirIndex = [dirData.isdir];
    fileList = {dirData(~dirIndex).name}';
    imageindexes = [];
    for i = 1:25
        load([num2str(i) '-saliencemap'],'fullmap');
        if any(any(isnan(fullmap)))
            disp('nans :(')
        end
        avgmap = avgmap+fullmap;
    end
end
imagesc(avgmap)
title('Average Salience Map')
axis off
%%
%---[3] Calculate Salience at Each Fixation---%
SMT_dir = 'C:\Users\seth.koenig\Documents\MATLAB\Scene Memory Task\'; %parent directory
load('allcortexfiles.mat'); %easily accessible data structure from above

imageX = 800; %horizontal pixel size of images
imageY = 600; %vertical pixel size of images

for crow = 12:size(allcortexfiles,1);
    cortexfile = allcortexfiles{crow,1};
    imset = allcortexfiles{crow,3};
    if imset < 10
        image_dir = [SMT_dir  'Image Sets\SMT00' num2str(imset) '\'];
    else
        image_dir = [SMT_dir  'Image Sets\SMT0' num2str(imset) '\'];
    end
    fixationSalience_and_significanceSMT(['C:\Users\seth.koenig\Documents\MATLAB\Scene Memory Task\'...
        'datafiles\' cortexfile '-fixation.mat'],image_dir,imageX,imageY);
end
%%
%---[4] Calculate Average Salience at Each Fixation Across Multiples Sets and Monkeys---%
SMT_dir = 'C:\Users\seth.koenig\Documents\MATLAB\Scene Memory Task\'; %parent directory

load('allcortexfiles.mat'); %easily accessible data structure from above

tags = cell(1,size(allcortexfiles,1)); %tags are unique to each monkey's names
for crow = 1:size(allcortexfiles,1);
    tags{crow} = allcortexfiles{crow,1}(1:2);
end
tags = unique(tags); %unique monkey identifiers

imageX = 800; %horizontal pixel size of images
imageY = 600; %vertical pixel size of images

minfix = cell(1,length(tags)); %buffer minimum number of fixations with larger number
for i = 1:length(minfix);
    minfix{i} = 100;
end

for crow = 1:size(allcortexfiles,1);
    cortexfile = allcortexfiles{crow,1};
    tag = find(ismember(tags,cortexfile(1:2)));
    load(['C:\Users\seth.koenig\Documents\MATLAB\Scene Memory Task\datafiles\'...
        cortexfile '-FixationStatistics.mat']);
    minfix{tag} = min([minfix{tag} size(statistics.numbervalues,2)]); %minimum number of good fixations
end

%cells by number of subjects,spacing,then lesion status (1 health, 2 lesion)
%spacing  0, 2, 4, 8, or 16 as indexes 1, 2, 3, 4, 5, respecitvely
novsal= cell(length(tags),5,2); %salience at fixations during novel presentations by spacing
famsal = cell(length(tags),5,2);%salience at fixations during familiar presentations by spacing
novimgI= cell(length(tags),5,2);%image intensity at fixations during novel presentations by spacing
famimgI = cell(length(tags),5,2);%image intensity at fixations during familiar presentations by spacing

%random is random so better estimate with more data
shuffsal= cell(length(tags),5,2); %salience at random fixations by spacing
shuffimgI= cell(length(tags),5,2);%image intensity at random fixations by spacing

count = ones(length(tags),5,2); %keep track of number of images per monkey, spacing, & lesion status

novnumberfixations = cell(length(tags),5,2); %novel presentation # of fixations. Will use for statistics testing later
famnumberfixations = cell (length(tags),5,2); %familiar presentation # of fixations. Will use for statistics testing later

%for loops combines the data
for crow = 1:size(allcortexfiles,1);
    cortexfile = allcortexfiles{crow,1};
    
    lesionstatus = allcortexfiles{crow,2};
    if strcmpi(lesionstatus,'H')%H for healthy control
        ls = 1;%if healhty place in z-column 1
    else
        ls = 2;%if 'L' for lesion place in z-column 1
    end
    
    tag = find(ismember(tags,cortexfile(1:2)));
    
    load([SMT_dir 'datafiles\' cortexfile '-FixationStatistics.mat']);%get salience and image intensity at fixations
    load([SMT_dir 'datafiles\' cortexfile '-fixation.mat'],'images'); %get image presentation order
    
    %grab salience values at fixations and reshape from 3D to 2D array
    saldata = shuffunshuffdata{2}{3}(:,1,:);
    saldata = reshape(saldata,[size(saldata,1),size(saldata,3)]);
    salshuffdata = shuffunshuffdata{1}{3}(:,1,:);
    salshuffdata  = reshape(salshuffdata ,[size(salshuffdata,1),size(salshuffdata,3)]);
    
    %grab image intesnity values at fixations and reshape from 3D to 2D array
    imgIdata = shuffunshuffdata{2}{3}(:,2,:);
    imgIdata = reshape(imgIdata,[size(imgIdata,1),size(imgIdata,3)]);
    imgIshuffdata = shuffunshuffdata{1}{3}(:,2,:);
    imgIshuffdata  = reshape(imgIshuffdata ,[size(imgIshuffdata,1),size(imgIshuffdata,3)]);
    
    for img = 1:25;
        ind = find(images == img);
        
        if length(ind) < 2; %use knowledge of spacings but ideally could fix, contructed in a structured manner
            if img <= 5
                spaceing = 5;
            elseif img <= 10
                spaceing = 4;
            elseif img <= 15
                spaceing = 3;
            elseif img <= 20
                spaceing = 2;
            else
                spaceing = 1;
            end
        else
            spaceing = nextpow2(diff(ind)-1)+1; %fancy way of chaning spacing as 0,2,4,8,16 into indexes 1,2,3,4,5
        end
        
        if ~isempty(ind)

            novsal{tag,spaceing,ls}(count(tag,spaceing,ls),:) = saldata(ind(1),1:minfix{tag});
            novimgI{tag,spaceing,ls}(count(tag,spaceing,ls),:) = imgIdata(ind(1),1:minfix{tag});
            
            if length(ind) == 2;%if novel and familiar presentations
                famsal{tag,spaceing,ls}(count(tag,spaceing,ls),:)= saldata(ind(2),1:minfix{tag});
                famimgI{tag,spaceing,ls}(count(tag,spaceing,ls),:)= imgIdata(ind(2),1:minfix{tag});
                
                %Combine shuffled data regardless of novel or familiar presentation
                cnt = 2*count(tag,spaceing,ls)-1;% double since 2x as many points plus odd index
                shuffsal{tag,spaceing,ls}(cnt,:) =  salshuffdata(ind(2),1:minfix{tag});
                shuffimgI{tag,spaceing,ls}(cnt,:) = imgIshuffdata(ind(2),1:minfix{tag});
                
                famnumberfixations{tag,spaceing,ls}(count(tag,spaceing,ls)) = sum(~isnan(saldata(ind(2),:)));
            else %if not want to preserve data structure so use NaNs
                famsal{tag,spaceing,ls}(count(tag,spaceing,ls),:)= NaN(1,minfix{tag});
                famimgI{tag,spaceing,ls}(count(tag,spaceing,ls),:)= NaN(1,minfix{tag});
                
                %Combine shuffled data regardless of novel or familiar presentation
                cnt = 2*count(tag,spaceing,ls)-1;% double since 2x as many points plus odd index
                shuffsal{tag,spaceing,ls}(cnt,:) =  NaN(1,minfix{tag});
                shuffimgI{tag,spaceing,ls}(cnt,:) = NaN(1,minfix{tag});
                
                famnumberfixations{tag,spaceing,ls}(count(tag,spaceing,ls)) = NaN;
            end
            
            cnt = 2*count(tag,spaceing,ls); %double since 2x as many points
            shuffsal{tag,spaceing,ls}(cnt,:) =  salshuffdata(ind(1),1:minfix{tag});
            shuffimgI{tag,spaceing,ls}(cnt,:) = imgIshuffdata(ind(1),1:minfix{tag});
            
            novnumberfixations{tag,spaceing,ls}(count(tag,spaceing,ls)) = sum(~isnan(saldata(ind(1),:)));
        else %perserve data structure fill with NaNs
            novsal{tag,spaceing,ls}(count(tag,spaceing,ls),:) = NaN(1,minfix{tag});
            novimgI{tag,spaceing,ls}(count(tag,spaceing,ls),:) = NaN(1,minfix{tag});
            famsal{tag,spaceing,ls}(count(tag,spaceing,ls),:)= NaN(1,minfix{tag});
            famimgI{tag,spaceing,ls}(count(tag,spaceing,ls),:)= NaN(1,minfix{tag});
            cnt = 2*count(tag,spaceing,ls)-1;% double since 2x as many points plus odd index
            shuffsal{tag,spaceing,ls}(cnt,:) =  NaN(1,minfix{tag});
            shuffimgI{tag,spaceing,ls}(cnt,:) = NaN(1,minfix{tag});
            cnt = 2*count(tag,spaceing,ls); %double since 2x as many points
            shuffsal{tag,spaceing,ls}(cnt,:) =  NaN(1,minfix{tag});
            shuffimgI{tag,spaceing,ls}(cnt,:) = NaN(1,minfix{tag});
            novnumberfixations{tag,spaceing,ls}(count(tag,spaceing,ls)) = NaN;
        end
        count(tag,spaceing,ls) = count(tag,spaceing,ls)+1;
    end
end
count = count-1; %to remove extra count

data_at_fixations.novel_salience = novsal;
data_at_fixations.familiar_salience = famsal;
data_at_fixations.novel_imageintensity = novimgI;
data_at_fixations.familiar_imageintensity = famimgI;
data_at_fixations.number_images = count;
data_at_fixations.shuffled_salinece = shuffsal;
data_at_fixations.shuffled_imageintensity = shuffimgI;
data_at_fixations.format = {'Row: monkey','Col: spacing','Z-Col 1: healthy','Z-Col 2: lesion'};


% now that the data is combined by monkey, spaceing, and lesion status we can plot and do statistical tests

%1st test salience at fixations during the novel presentaion is greater than chance for all fixations across monkeys

%group novel salience data together by monkey across spacing
minnumfix = min(cell2mat(minfix));%minimum number of good fixations from all monkeys
mediannumfixs = floor(nanmedian(nanmedian(nanmedian(cellfun(@nanmedian,novnumberfixations)))));
novsaltag = cell(1,length(tags));
shuffsaltag = cell(1,length(tags));
for tag = 1:size(novsal,1);
    for space = 1:size(novsal,2);
        for ls = 1:2;
            if ~isempty(novsal{tag,space,ls});
                novsaltag{tag} = [novsaltag{tag}; novsal{tag,space,ls}];
                shuffsaltag{tag} = [shuffsaltag{tag}; shuffsal{tag,space,ls}];
            end
        end
    end
end

%group novel salience data together across all monkeys and spaceings of images
novsalall = [];
shuffsalall = [];
for tag = 1:size(novsaltag,2);
    novsalall = [novsalall; novsaltag{tag}(:,1:minnumfix)];
    shuffsalall = [shuffsalall; shuffsaltag{tag}(:,1:minnumfix)];
end

% Significance test for salience at fixations during novel presentation
% by fixation number across all monkeys and spaceings of images
allnovsalpvalues = NaN(1,size(novsalall,2));
allnovsalcI = NaN(1,size(novsalall,2));
shuffleddata = shuffsalall(1:end);
shuffleddata(isnan(shuffleddata)) = [];
for fix = 1:size(novsalall,2);
    [~,p,ci] = ztest(shuffleddata,nanmean(novsalall(:,fix)),std(shuffleddata),0.05,'left');
    allnovsalpvalues(fix) = p;
    allnovsalcI(fix) = ci(2);
end

% Significance test for salience at fixations during novel presentation
% by fixation number and by money monkeys across all spaceings of images
novsalpvalues = NaN(size(novsaltag,2),mediannumfixs);
novsalcI = NaN(size(novsaltag,2),mediannumfixs);
for tag = 1:size(novsaltag,2);
    for fix = 1:mediannumfixs
        [~,p,ci] = ztest(shuffleddata,nanmean(novsaltag{tag}(:,fix)),std(shuffleddata),0.05,'left');
        novsalpvalues(tag,fix) = p;
        novsalcI (tag,fix) = ci(2);
    end
end

figure
hold all
plot(nanmean(novsalall),'linewidth',5);
for tag = 1:size(novsaltag,2);
    plot(nanmean(novsaltag{tag}))
end
plot(1:mediannumfixs,allnovsalcI(1:mediannumfixs),'k--')
hold off
xlabel('Fixation Number')
ylabel('Salience')
box off
title('Salience at Fixation Locations during Novel Presentation across All Monkeys')
xlim([0 mediannumfixs+1]);
legend(['All monkeys',tags,'Chance Salience'])

stats_at_novel_fixations.salience_cI = novsalcI;
stats_at_novel_fixations.salience_pvalues = novsalpvalues;
stats_at_novel_fixations.combined_salinece_cI = allnovsalcI;
stats_at_novel_fixations.combined_salinece_cI = allnovsalpvalues;
stats_at_novel_fixations.format = {'Row: Monkey','Col: fixation #'};
stats_at_novel_fixations.combined_format = {'Column: fixation #'};

% 2nd test if salience at fixations during the novel presentaion is greater than the
% salinece at fixations during the familiar presentation. Will test on individual basis,
% spacing, and lesion status

averagesalience = cell(2,size(novsal,2));
novel_vs_familiar_pvalues = NaN(length(tags),size(novsal,2),2); %monkey by spaceing
for tag = 1:size(novsal,1);
    for ls = 1:2;
        if ~isempty(novsal{tag,1,ls})
            figure
            for space = 1:size(novsal,2);
                subplot(3,2,space)
                hold on
                plot(nanmean(novsal{tag,space,ls}(:,1:mediannumfixs)),'b')
                plot(nanmean(famsal{tag,space,ls}(:,1:mediannumfixs)),'r');
                xlabel('Fixation Number')
                ylabel('Salience')
                if space == 1;
                    title('Presented 0 images apart')
                else
                    title(['Presented ' num2str(2^(space-1)) ' images apart'])
                end
                hold off
                averagesalience{1,space} = [averagesalience{1,space}; novsal{tag,space,ls}(:,1:mediannumfixs)];
                averagesalience{2,space} = [averagesalience{2,space}; famsal{tag,space,ls}(:,1:mediannumfixs)];
            end
            
            subplot(3,2,6)
            hold on
            for space = 1:size(novsal,2);
                [~,p] = ttest2(averagesalience{1,space}(2:10),averagesalience{2,space}(2:10));
                novel_vs_familiar_pvalues(tag,space,ls) = p;
                numpoints = sqrt(sum(~isnan(averagesalience{1,space}(2:10))));
                errorbar(2*space-1,nanmean(averagesalience{1,space}(2:10)),...
                    nanstd(averagesalience{1,space}(2:10))/numpoints,'b')
                errorbar(2*space,nanmean(averagesalience{2,space}(2:10)),...
                    nanstd(averagesalience{2,space}(2:10))/numpoints,'r')
                if p < 0.05
                    plot(2*space-0.5,0.5,'k*')
                end
            end
            hold off
            set(gca,'XTick',1.5:2:11.5)
            set(gca,'XTickLabel',{'0','2','4','8','16'})
            xlabel('Number of images in between novel and familiar presentations')
            ylabel('Salience')
            title('Salience across all 2-10th fixations')
            legend('Novel','Familiar')
            
            if ls == 1;
                subtitle(['Salience at Fixations during Novel and Familiar Presentations for ' ...
                    tags{tag} ' as control']);
            else
                subtitle(['Salience at Fixations during Novel and Familiar Presentations for ' ...
                    tags{tag} ' with lesion']);
            end
        end
    end
end

%calculate if salience at fixations during familiar presentations is significantly greater than chance
familiar_salience_pvalues = NaN(length(tags),size(novsal,2),2); %monkey by spaceing, by lesion status
familiar_salience_CI = NaN(length(tags),size(novsal,2),2); %monkey by spaceing, by lesion status
familiar_salience_means = cell(length(tags),size(novsal,2),2); %monkey by spaceing, by lesion status
familiar_salience_std = cell(length(tags),size(novsal,2),2); %monkey by spaceing, by lesion status
for tag = 1:size(famsal,1);
    for space = 1:size(famsal,2);
        for ls = 1:size(famsal,3);
            if ~isempty(famsal{tag,space,ls})
                for fix = 1:mediannumfixs
                    [~,p,ci] = ztest(shuffleddata,nanmean(famsal{tag,space,ls}(:,fix)),std(shuffleddata),0.05,'left');
                    familiar_salience_pvalues(tag,space,ls) = p;
                    familiar_salience_CI(tag,space,ls)= ci(2);
                end
                familiar_salience_means{tag,space,ls} = nanmean(famsal{tag,space,ls}(1:mediannumfixs));
                familiar_salience_std{tag,space,ls} = nanstd(famsal{tag,space,ls}(1:mediannumfixs));
            end
        end
    end
end

stats_at_familiar_fixations.pvalues = familiar_salience_means;
stats_at_familiar_fixations.CI = familiar_salience_CI;
stats_at_familiar_fixations.means = familiar_salience_means;
stats_at_familiar_fixations.stds = familiar_salience_std;
stats_at_familiar_fixations.format = {'Row: Monkey','Col: image spaceing',...
    'Z-Col 1: healthy','Z-Col 2: lesion'};


%Statistical significance on salience at novel vs familiar locations for pre vs post lesion
pre_postlesion_novel_vs_familiar_pvalues = NaN(length(tags),size(novsal,2),2); %monkey by spaceing by novel vs familiar
for tag = 1:size(novsal,1);
    if ~isempty(novsal{tag,1,1}) && ~isempty(novsal{tag,1,2})
        for space = 1:size(novsal,2)
            [~,p] = ttest2(novsal{tag,space,1},novsal{tag,space,1});%is there a difference during novel presentations?
            pre_postlesion_novel_vs_familiar_pvalues(tag,space,1) = p;
            [~,p] = ttest2(novsal{tag,space,1},novsal{tag,space,2});%is there a difference during familiar presentations?
            pre_postlesion_novel_vs_familiar_pvalues(tag,space,2) = p;
        end
        disp('Need to put more analysis here (line 536), a plot, etc.!!!')
    end
end

%Statistical analysis comparing heatlhy monkeys to lesion monkeys for salience at
% fixations during novel and familiar presentations
combinedcontrol = cell(2,size(novsal,2));
combinedlesion = cell(2,size(novsal,2));
for tag = 1:size(novsal,1);
    for space = 1:size(novsal,2);
        if ~isempty(novsal{tag,space,1})
            combinedcontrol{1,space} = [combinedcontrol{1,space};  novsal{tag,space,1}(:,1:mediannumfixs)];
            combinedcontrol{2,space} = [combinedcontrol{2,space};  famsal{tag,space,1}(:,1:mediannumfixs)];
        end
        if ~isempty(novsal{tag,space,2})
            combinedlesion{1,space} = [combinedlesion{1,space};  novsal{tag,space,2}(:,1:mediannumfixs)];
            combinedlesion{2,space} = [combinedlesion{2,space};  famsal{tag,space,2}(:,1:mediannumfixs)];
        end
    end
end

healthy_vs_lesion_pvalues = NaN(size(novsal,2),2); %row spaceing, column novelty: novel 1 then 2 familiar
healthy_vs_lesion_means = NaN(size(novsal,2),2,2); %1 for healthy 2 for lesion
healthy_vs_lesion_std = NaN(size(novsal,2),2,2);
healthy_vs_lesion_dmeans = NaN(size(novsal,2),2); %difference in mean values
healthy_vs_lesion_dstd = NaN(size(novsal,2),2);%difference std values
healthy_vs_lesion_dnumpoints = NaN(size(novsal,2),2);%difference std values
healthy_vs_lesion_numpoints = NaN(size(novsal,2),2,2);
for space = 1:size(combinedcontrol,2);
    for novelty = 1:size(combinedcontrol,1) %novel vs familiar
        [~,p] = ttest2(combinedcontrol{novelty,space}(1:end),combinedlesion{novelty,space}(1:end));
        healthy_vs_lesion_pvalues(space,novelty) = p;
        
        %healthy monkeys
        healthy_vs_lesion_means(space,novelty,1) = nanmean(combinedcontrol{novelty,space}(1:end));
        healthy_vs_lesion_std(space,novelty,1) = nanstd(combinedcontrol{novelty,space}(1:end));
        healthy_vs_lesion_numpoints(space,novelty,1) = sum(~isnan(combinedcontrol{novelty,space}(1:end)));
        
        %lesion monkeys
        healthy_vs_lesion_means(space,novelty,2) = nanmedian(combinedlesion{novelty,space}(1:end));
        healthy_vs_lesion_std(space,novelty,2) = nanstd(combinedlesion{novelty,space}(1:end));
        healthy_vs_lesion_numpoints(space,novelty,2) = sum(~isnan(combinedlesion{novelty,space}(1:end)));
    end
    %Change in salience from novel to familiar presentation...as control
    healthy_vs_lesion_dmeans(space,1) = nanmean(combinedcontrol{1,space}(1:end)-...%if no repeat then NaN
        combinedcontrol{2,space}(1:end));
    healthy_vs_lesion_dstd(space,1) = nanstd(combinedcontrol{1,space}(1:end)-...%if no repeat then NaN
        combinedcontrol{2,space}(1:end));
    healthy_vs_lesion_dnumpoints(space,1) = sum(~isnan(combinedcontrol{1,space}(1:end)-...%if no repeat then NaN
        combinedcontrol{2,space}(1:end)));
    healthy_vs_lesion_dmeans(space,2) = nanmean(combinedlesion{1,space}(1:end)-...%if no repeat then NaN
        combinedlesion{2,space}(1:end));
    healthy_vs_lesion_dstd(space,2) = nanstd(combinedlesion{1,space}(1:end)-...%if no repeat then NaN
        combinedlesion{2,space}(1:end));
    healthy_vs_lesion_dnumpoints(space,2) = sum(~isnan(combinedlesion{1,space}(1:end)-...%if no repeat then NaN
        combinedlesion{2,space}(1:end)));
end

healthy_vs_lesion.means = healthy_vs_lesion_means;
healthy_vs_lesion.stds = healthy_vs_lesion_std;
healthy_vs_lesion.numpoints = healthy_vs_lesion_numpoints;
healthy_vs_lesion.pvalues = healthy_vs_lesion_pvalues;
healthy_vs_lesion.form = {'Row: image spaceing','Col 1 novel',...
    'Col 2: familiar','Z-Col 1 healthy','Z-Col 2: lesion'};


figure
subplot(1,2,1)
hold on
errorbar(healthy_vs_lesion_means(:,1,1),...
    healthy_vs_lesion_std(:,1,1)./sqrt(healthy_vs_lesion_numpoints(:,1,1)),'b');
errorbar(healthy_vs_lesion_means(:,1,2),...
    healthy_vs_lesion_std(:,1,2)./sqrt(healthy_vs_lesion_numpoints(:,1,2)),'r');
hold off
legend('Healthy Novel','Lesion Novel')
set(gca,'XTick',1:5)
set(gca,'XTickLabel',{'0','2','4','8','16'})
xlabel('Number of images in between novel and familiar presentations')
ylabel('Salience')
ylim([0.30 0.44])
title('Salience at Fixations during Novel Presetions: Control vs. Lesion Monkeys')

subplot(1,2,2)
hold on
errorbar(healthy_vs_lesion_means(:,2,1),...
    healthy_vs_lesion_std(:,2,1)./sqrt(healthy_vs_lesion_numpoints(:,2,1)),'b');
errorbar(healthy_vs_lesion_means(:,2,2),...
    healthy_vs_lesion_std(:,2,2)./sqrt(healthy_vs_lesion_numpoints(:,2,2)),'r');
hold off
legend('Healthy Familiar','Lesion Familiar')
set(gca,'XTick',1:5)
set(gca,'XTickLabel',{'0','2','4','8','16'})
xlabel('Number of images in between novel and familiar presentations')
ylabel('Salience')
ylim([0.30 0.44])
title('Salience at Fixations during Familiar Presetions: Control vs. Lesion Monkeys')

figure
hold on
errorbar(healthy_vs_lesion_dmeans(:,1),...
    healthy_vs_lesion_dstd(:,1)./sqrt(healthy_vs_lesion_dnumpoints(:,1)),'b');
errorbar(healthy_vs_lesion_dmeans(:,2),...
    healthy_vs_lesion_dstd(:,2)./sqrt(healthy_vs_lesion_dnumpoints(:,2)),'r');
hold off
legend('Healthy','Lesion','location','NorthEastOutside')
set(gca,'XTick',1:5)
set(gca,'XTickLabel',{'0','2','4','8','16'})
xlabel('Number of images in between novel and familiar presentations')
ylabel('\Delta Salience')
title('Change in Salience at Fixations from Novel to Familiar Presentations')


variables = {
    'All structure variables have a format field that tells you how to read the data';
    '';
    'data_at_fixations: salience and and image intensity values at fixation locations';
    'as well as shuffled (random) locations';
    '';
    'stats_at_novel_fixations: salience pvalues, means, stds at fixations during novel presentions';
    'stats_at_familiar_fixations: salience pvalues, means, stds at fixations during familiar presentions';
    '';
    'healthy_vs_lesion: salience pvalues, means, stds at fixations comparing healthy to lesion monkeys';
    '';
    'pre_postlesion_novel_vs_familiar_pvalues: pvalues testing difference in salience at fixations'
    'during novel vs familiar presentations if a monkey had pre and post lesion data';
    '';
    'minfix: minimum number of valid fixations by monkey';
    '';
    'mediannumfix: median number of fixations for stastical and plotting purposes';
    '';
    'tags: monkey initials';
    ''
    'novel_vs_familiar_pvalues: statstical analysis results testing if there was a diffrence';
    'in salience values at fixations during novel vs familiar presentations. Formatted';
    'with row corresponding to monkey and row corresponding to image spaceing';
    };

save([SMT_dir 'datafiles\Combined_Salience_Data.mat'],'data_at_fixations',...
    'stats_at_novel_fixations','stats_at_familiar_fixations','healthy_vs_lesion',...
    'pre_postlesion_novel_vs_familiar_pvalues','minfix','mediannumfixs','tags',...
    'novel_vs_familiar_pvalues','variables')
%%
%---[5] Calculate Similarity in Fixation Locations with KL Divergence---%
SMT_dir = 'C:\Users\seth.koenig\Documents\MATLAB\Scene Memory Task\'; %parent directory

load('allcortexfiles.mat'); %easily accessible data structure from above

tags = cell(1,size(allcortexfiles,1)); %tags are unique to each monkey's names
for crow = 1:size(allcortexfiles,1);
    tags{crow} = allcortexfiles{crow,1}(1:2);
end
tags = unique(tags); %unique monkey identifiers

f = fspecial('gaussian',[256,256],24); %1 dva 2D gaussian filter
binsize = 25;
imageX = 800; %horizontal pixel size of images
imageY = 600; %vertical pixel size of images

overlap = cell(length(tags),5,2);
shuffledoverlap = cell(length(tags),5,2);
spaceingind = ones(length(tags),2,6,5);
for crow = 1:size(allcortexfiles,1);
    cortexfile = allcortexfiles{crow,1};
    tag = find(ismember(tags,cortexfile(1:2)));
    lesionstatus = allcortexfiles{crow,2};
    if strcmpi(lesionstatus,'H')%H for healthy control
        ls = 1;%if healhty place in z-column 1
    else
        ls = 2;%if 'L' for lesion place in z-column 1
    end
    
    load([SMT_dir 'datafiles\' cortexfile '-fixation.mat'])
    
    for img = 1:25;
        ind = find(images == img);
        if length(ind) == 2
            spaceing = nextpow2(diff(ind)-1)+1;
            
            novelmap = cell(1,6);
            familiarmap = cell(1,6);
            shuffnovelmap = cell(1,6);
            shufffamiliarmap = cell(1,6);
            for i = 1:6;
                novelmap{i} = zeros(imageY,imageX);
                familiarmap{i} = zeros(imageY,imageX);
                shuffnovelmap{i} = zeros(imageY,imageX);
                shufffamiliarmap{i} = zeros(imageY,imageX);
            end
            
            %break fixations into groups of 5 so we can get order and location information
            novelfixations = fixationstats{ind(1)}.fixations;
            familiarfixations = fixationstats{ind(2)}.fixations;
            
            if ~isempty(novelfixations)
                if novelfixations(1,1) > imageX/2-100 && novelfixations(1,1) < imageX/2+100 &&...
                        novelfixations(2,1) < imageY/2+100 && novelfixations(2,1) > imageY/2-100
                    novelfixations(:,1) = [];
                    novelfixations(:,1) = [];
                end
            end
            
            if ~isempty(familiarfixations)
                if familiarfixations(1,1) > imageX/2-100 && familiarfixations(1,1) < imageX/2+100 &&...
                        familiarfixations(2,1) < imageY/2+100 && familiarfixations(2,1) > imageY/2-100
                    familiarfixations(:,1) = [];
                    familiarfixations(:,1) = [];
                end
            end
            
            maxfixations = min(size(novelfixations,2),size(familiarfixations,2));
            maxfiations(maxfixations > 25) = 25;
            maxfixations = maxfixations-rem(maxfixations,5);
            for fix = 1:maxfixations;
                fixx = round(novelfixations(1,fix));
                fixx(fixx < 1) = 1; fixx(fixx > imageX) = imageX;
                fixy = imageY - round(novelfixations(2,fix));
                fixy(fixy < 1) = 1; fixy(fixy > imageY) = imageY;
                randx = randi(imageX);randx(randx < 1) = 1; %can occassionaly get 0 so have to correct
                randy = randi(imageY);randy(randy < 1) = 1;%can occassionaly get 0 so have to correct
                if fix <= 5
                    novelmap{1}(fixy,fixx) = novelmap{1}(fixy,fixx) + 1;
                    shuffnovelmap{1}(randy,randx) = shuffnovelmap{1}(randy,randx) + 1;
                elseif fix <= 10
                    novelmap{2}(fixy,fixx) = novelmap{2}(fixy,fixx) + 1;
                    shuffnovelmap{2}(randy,randx) = shuffnovelmap{2}(randy,randx) + 1;
                elseif fix <= 15
                    novelmap{3}(fixy,fixx) = novelmap{3}(fixy,fixx) + 1;
                    shuffnovelmap{3}(randy,randx) = shuffnovelmap{3}(randy,randx) + 1;
                elseif fix <= 20
                    novelmap{4}(fixy,fixx) = novelmap{4}(fixy,fixx) + 1;
                    shuffnovelmap{4}(randy,randx) = shuffnovelmap{4}(randy,randx) + 1;
                elseif fix <= 25 && (length(novelfixations) >= 25)
                    novelmap{5}(fixy,fixx) = novelmap{5}(fixy,fixx) + 1;
                    shuffnovelmap{5}(randy,randx) = shuffnovelmap{5}(randy,randx) + 1;
                end
            end
            
            for fix = 1:maxfixations;
                fixx = round(familiarfixations(1,fix));
                fixx(fixx < 1) = 1; fixx(fixx > imageX) = imageX;
                fixy = imageY - round(familiarfixations(2,fix));
                fixy(fixy < 1) = 1; fixy(fixy > imageY) = imageY;
                randx = randi(imageX);
                randx(randx < 1) = 1; %can occassionaly get 0 so have to correct
                randy = randi(imageY);
                randy(randy < 1) = 1; %can occassionaly get 0 so have to correct
                if fix <= 5
                    familiarmap{1}(fixy,fixx) = familiarmap{1}(fixy,fixx) + 1;
                    shufffamiliarmap{1}(randy,randx) = shufffamiliarmap{1}(randy,randx) + 1;
                elseif fix <= 10
                    familiarmap{2}(fixy,fixx) = familiarmap{2}(fixy,fixx) + 1;
                    shufffamiliarmap{2}(randy,randx) = shufffamiliarmap{2}(randy,randx) + 1;
                elseif fix <= 15
                    familiarmap{3}(fixy,fixx) = familiarmap{3}(fixy,fixx) + 1;
                    shufffamiliarmap{3}(randy,randx) = shufffamiliarmap{3}(randy,randx) + 1;
                elseif fix <= 20
                    familiarmap{4}(fixy,fixx) = familiarmap{4}(fixy,fixx) + 1;
                    shufffamiliarmap{4}(randy,randx) = shufffamiliarmap{4}(randy,randx) + 1;
                elseif fix <= 25 && (length(novelfixations) >= 25)
                    familiarmap{5}(fixy,fixx) = familiarmap{5}(fixy,fixx) + 1;
                    shufffamiliarmap{5}(randy,randx) = shufffamiliarmap{5}(randy,randx) + 1;
                end
            end
            
            % combine fixations 1-25 into a single map so we get just location overlap
            for i = 1:5
                novelmap{6} = novelmap{6}+novelmap{i};
                shuffnovelmap{6} = shuffnovelmap{6}+shuffnovelmap{i};
                familiarmap{6} = familiarmap{6}+familiarmap{i};
                shufffamiliarmap{6} = shufffamiliarmap{6}+shufffamiliarmap{i};
            end
            
            % For loop below...
            % 1. convoles fixation matrix with 1 dva (24 pixel) gaussian filter
            % to account for variability in fixation locations and eye tracking error
            % 2. Bins convolved matrix
            % 3. removes 0's and replaces with minimum defined value in matlab eps (2 ^-52)
            % 4. creates fixation pdf by dividing matrix by the sum of the matrix
            for i = 1:6;
                novelmap{i} = imfilter(novelmap{i},f);
                shuffnovelmap{i} = imfilter(shuffnovelmap{i},f);
                novelmap{i} = bin2(novelmap{i},binsize,binsize);
                shuffnovelmap{i} = bin2(shuffnovelmap{i},binsize,binsize);
                novelmap{i}(novelmap{i} == 0) = eps;
                shuffnovelmap{i}(shuffnovelmap{i} == 0) = eps;
                novelmap{i} = novelmap{i}./sum(sum(novelmap{i}));
                shuffnovelmap{i} = shuffnovelmap{i}./sum(sum(shuffnovelmap{i}));
                
                familiarmap{i} = imfilter(familiarmap{i},f);
                shufffamiliarmap{i} = imfilter(shufffamiliarmap{i},f);
                familiarmap{i} = bin2(familiarmap{i},binsize,binsize);
                shufffamiliarmap{i} = bin2(shufffamiliarmap{i},binsize,binsize);
                familiarmap{i}(familiarmap{i} == 0) = eps;
                shufffamiliarmap{i}(shufffamiliarmap{i} == 0) = eps;
                familiarmap{i} = familiarmap{i}./sum(sum(familiarmap{i}));
                shufffamiliarmap{i} = shufffamiliarmap{i}./sum(sum(shufffamiliarmap{i}));
            end
            
            %calculate overlap using KL divergence. Close more similar PDFs have smaller KL divergence (distance in bits)
            for i = 1:size(novelmap,2)
                if ~all(all(novelmap{i} == mean(mean(novelmap{i})))) %if maps are not empty
                    overlap{tag,spaceing,ls}{i}(spaceingind(tag,ls,i,spaceing)) = ...
                        sum(sum(log2(novelmap{i}./familiarmap{i}).*novelmap{i})) + ...
                        sum(sum(log2(familiarmap{i}./novelmap{i}).*familiarmap{i}));
                    
                    shuffledoverlap{tag,spaceing,ls}{i}(spaceingind(tag,ls,i,spaceing)) = ...
                        sum(sum(log2(shuffnovelmap{i}./shufffamiliarmap{i}).*shuffnovelmap{i})) + ...
                        sum(sum(log2(shufffamiliarmap{i}./shuffnovelmap{i}).*shufffamiliarmap{i}));
                else
                    overlap{tag,spaceing,ls}{i}(spaceingind(tag,ls,i,spaceing)) = NaN;
                    shuffledoverlap{tag,spaceing,ls}{i}(spaceingind(tag,ls,i,spaceing)) = NaN;
                end
                spaceingind(tag,ls,i,spaceing) = spaceingind(tag,ls,i,spaceing)+1;
            end
        else %want to preserve data structure so fill with NaNs
            if img <= 5
                spaceing = 5;
            elseif img <= 10
                spaceing = 4;
            elseif img <= 15
                spaceing = 3;
            elseif img <= 20
                spaceing = 2;
            else
                spaceing = 1;
            end
            for i = 1:6
                overlap{tag,spaceing,ls}{i}(spaceingind(tag,ls,i,spaceing)) = NaN;
                shuffledoverlap{tag,spaceing,ls}{i}(spaceingind(tag,ls,i,spaceing)) = NaN;
                spaceingind(tag,ls,i,spaceing) = spaceingind(tag,ls,i,spaceing)+1;
            end
        end
    end
end

combinedshuffled = cell(1,6);
for ls = 1:size(overlap,3);
    for tag = 1:size(overlap,1);
        for space = 1:size(overlap,2);
            for fixgroup = 1:6; % for groups of fixations 1-5,6-10,11-15,16-20,21-25,all/1-25
                if ~isempty(shuffledoverlap{tag,space,ls})
                    combinedshuffled{fixgroup} = [combinedshuffled{fixgroup} shuffledoverlap{tag,space,ls}{fixgroup}];
                end
            end
        end
    end
end

meanKLdivergence = NaN(size(overlap,1),size(overlap,2),6,2);
stdKLdivergence = NaN(size(overlap,1),size(overlap,2),6,2);
numpointsKLdivergence = NaN(size(overlap,1),size(overlap,2),6,2);
pvaluesKLdivergence = NaN(size(overlap,1),size(overlap,2),6,2); %differnce from chance
for ls = 1:size(overlap,3);
    for tag = 1:size(overlap,1);
        for space = 1:size(overlap,2);
            if ~isempty(overlap{tag,space,ls});
                for fixgroup = 1:6; % for groups of fixations 1-5,6-10,11-15,16-20,21-25,all/1-25
                    if ~all(isnan(overlap{tag,space,ls}{fixgroup}))
                        [~,p] = kstest2(overlap{tag,space,ls}{fixgroup},combinedshuffled{fixgroup});
                        pvaluesKLdivergence(tag,space,fixgroup) = p;
                        meanKLdivergence(tag,space,fixgroup,ls) = nanmean(overlap{tag,space,ls}{fixgroup});
                        stdKLdivergence(tag,space,fixgroup,ls) = nanstd(overlap{tag,space,ls}{fixgroup});
                        numpointsKLdivergence(tag,space,fixgroup,ls) = sum(~isnan(overlap{tag,space,ls}{fixgroup}));
                    end
                end
            end
        end
    end
end

for tag = 1:size(meanKLdivergence,1)
    figure
    hold on
    for ls = 1:size(meanKLdivergence,4);
        for space = 1:size(meanKLdivergence,2)
            subplot(3,2,space)
            
            mn = reshape(meanKLdivergence(tag,space,:,ls),[1,6]);
            sem = reshape(stdKLdivergence(tag,space,:,ls),[1,6])./reshape(numpointsKLdivergence(tag,space,:,ls),[1,6]);
            
            hold on
            errorbar(mn,sem,'-')
            plot(0.5:5.5,ones(1,6).*nanmean(cellfun(@nanmean,combinedshuffled)),'k--')
            plot(5.5:6.5,ones(1,2).*nanmean(combinedshuffled{6}),'k--')
            hold off
            
            set(gca,'XTick',[1:6])
            set(gca,'XTickLabel',{'1-5','6-10','11-15','16-20','21-25','1-25'});
            xlabel('Group of fixations')
            ylabel('KL Divergence (Distance in bits)')
            title([num2str(num2str(2^(space-1))) ' between Novel and Familiar Images Presentations'])
        end
        if ls == 1
            subtitle(['Overlap in Fixations Locations during Novel and Familiar Presentations for ' tags{tag} 'with lesion']);
        else
            subtitle(['Overlap in Fixations Locations during Novel and Familiar Presentations for ' tags{tag} 'as healthy control']);
        end
    end
end

legendlabels = {};
figure
hold all
for tag = 1:size(meanKLdivergence,1)
    for ls = 1:size(meanKLdivergence,4);
        if all(~isnan(meanKLdivergence(tag,:,6,ls)));
            mn = reshape(meanKLdivergence(tag,:,6,ls),[1,5]);
            sem = reshape(stdKLdivergence(tag,:,6,ls),[1,5])./reshape(numpointsKLdivergence(tag,:,6,ls),[1,5]);
            errorbar(mn,sem,'-')
            if ls == 1
                legendlabels = [legendlabels {[tags{tag} ' as healthy control']}];
            else
                legendlabels = [legendlabels {[tags{tag} ' w/ lesion']}];
            end
        end
    end
end
set(gca,'XTick',[1:5])
set(gca,'XTickLabel',{'0','2','4','8','16'});
xlabel('Number of Images in between Novel and Familiar Presentations')
ylabel('KL Divergence (Distance in bits)')
title('Overlap in Fixations Locations during Novel and Familiar Presentations')
legend(legendlabels)

healthyKL = cell(6,5);
lesionKL = cell(6,5);
for ls = 1:size(overlap,3);
    for tag = 1:size(overlap,1);
        for space = 1:size(overlap,2);
            if ~isempty(overlap{tag,space,ls});
                for fixgroup = 1:6; % just for all fixations 1-25
                    if ls == 1;
                        healthyKL{fixgroup,space} = [healthyKL{fixgroup,space} overlap{tag,space,ls}{fixgroup}];
                    else
                        lesionKL{fixgroup,space} = [lesionKL{fixgroup,space} overlap{tag,space,ls}{fixgroup}];
                    end
                end
            end
        end
    end
end

fixation_groups = {'1-5','6-10','11-15','16-20','21-25','1-25'};
figure
for map = 1:6
    subplot(2,3,map)
    hold on
    plot(cellfun(@nanmean,healthyKL(map,:)),'b')
    plot(cellfun(@nanmean,lesionKL(map,:)),'r')
    hold off
    legend('Heatlhy Control','Lesion')
    set(gca,'XTick',[1:5])
    set(gca,'XTickLabel',{'0','2','4','8','16'});
    xlabel('Number of Images in between Novel and Familiar Presentations')
    ylabel('KL Divergence (Distance in bits)')
    title(['Fixations ' fixation_groups{map}])
end
subtitle('Overlap in Fixations Locations during Novel and Familiar Presentations-healthy vs lesion')
%%
%---[6] Calculate Viewing Behavior Statistics (e.g. fixation duration and saccade amplitude)---%
SMT_dir = 'C:\Users\seth.koenig\Documents\MATLAB\Scene Memory Task\'; %parent directory

load('allcortexfiles.mat'); %easily accessible data structure from above

imageX = 800; %horizontal pixel size of images
imageY = 600; %vertical pixel size of images
SAMPRATE = 5;
for crow = 1:size(allcortexfiles,1);
    cortexfile = allcortexfiles{crow,1};
    getViewingBehaviorSMT([SMT_dir 'datafiles\' cortexfile '-fixation.mat'],SAMPRATE,imageX,imageY)
end
%%
%--- [7] Combine Viewing Behavior Across Monkeys and Lesion Status---%
SMT_dir = 'C:\Users\seth.koenig\Documents\MATLAB\Scene Memory Task\'; %parent directory

load('allcortexfiles.mat'); %easily accessible data structure from above

tags = cell(1,size(allcortexfiles,1)); %tags are unique to each monkey's names
for crow = 1:size(allcortexfiles,1);
    tags{crow} = allcortexfiles{crow,1}(1:2);
end
tags = unique(tags); %unique monkey identifiers

f = fspecial('gaussian',[256,256],24); %1 dva 2D gaussian filter
samprate = 5;

imageX = 800; %horizontal pixel size of images
imageY = 600; %vertical pixel size of images

allview = cell(length(tags),5,2);%monkeys by spacing by lesion status
for t = 1:length(tags)
    for space = 1:5
        for ls = 1:2
            %for novel image presentations
            allview{t,space,ls}.nov_densitymap = zeros(imageY,imageX); %PDF of fixation locations
            allview{t,space,ls}.nov_anglebtwfix = []; %angle between fixation locations
            allview{t,space,ls}.nov_sacangle_2fix = []; %saccade angle entering a fixation
            allview{t,space,ls}.nov_distbtwnfix = [];%distance between fixations
            allview{t,space,ls}.nov_fixduration = [];% fixaiton durations
            allview{t,space,ls}.nov_sacangle = [];%saccade angle leaving a fixation
            allview{t,space,ls}.nov_sacdist = [];%saccade arc length
            allview{t,space,ls}.nov_sacamplitude = [];%saccade amplitude
            allview{t,space,ls}.nov_sacduration = [];%saccade duration
            allview{t,space,ls}.nov_timebtwfix = [];%time between fixations
            
            %for familiar image presentations
            allview{t,space,ls}.fam_densitymap = zeros(imageY,imageX);
            allview{t,space,ls}.fam_anglebtwfix = [];
            allview{t,space,ls}.fam_sacangle_2fix = [];
            allview{t,space,ls}.fam_distbtwnfix = [];
            allview{t,space,ls}.fam_fixduration = [];
            allview{t,space,ls}.fam_sacangle = [];
            allview{t,space,ls}.fam_sacdist = [];
            allview{t,space,ls}.fam_sacamplitude = [];
            allview{t,space,ls}.fam_sacduration = [];
            allview{t,space,ls}.fam_timebtwfix = [];
        end
    end
end

spaceingind = ones(length(tags),5,2);
for crow = 1:size(allcortexfiles,1);
    cortexfile = allcortexfiles{crow,1};
    tag = find(ismember(tags,cortexfile(1:2)));
    lesionstatus = allcortexfiles{crow,2};
    if strcmpi(lesionstatus,'H')%H for healthy control
        ls = 1;%if healhty place in z-column 1
    else
        ls = 2;%if 'L' for lesion place in z-column 1
    end
    
    load([SMT_dir 'datafiles\' cortexfile '-fixation.mat'])
    load([SMT_dir 'datafiles\' cortexfile '-ViewingBehavior.mat'])
    
    for img = 1:25;
        ind = find(images == img);
        if length(ind) == 2
            spaceing = nextpow2(diff(ind)-1)+1;
            
            fixations = fixationstats{ind(1)}.fixations;
            for f = 2:size(fixations,2);  %ignore first since it's a fixaiton on the cross hair
                xx = round(fixations(1,f));
                yy = round(imageY-fixations(2,f));
                xx(xx < 1) = 1;
                yy(yy < 1)= 1;
                xx(xx > imageX) = imageX;
                yy(yy > imageY) = imageY;
                allview{tag,spaceing,ls}.nov_densitymap(yy,xx) = ...
                    allview{tag,spaceing,ls}.nov_densitymap(yy,xx)+1;
            end
            allview{tag,spaceing,ls}.nov_anglebtwfix = [allview{tag,spaceing,ls}.nov_anglebtwfix; ...
                anglebtwfix(ind(1),:)];
            allview{tag,spaceing,ls}.nov_sacangle_2fix = [allview{tag,spaceing,ls}.nov_sacangle_2fix; ...
                sacangle_2fix(ind(1),:)];
            allview{tag,spaceing,ls}.nov_distbtwnfix = [allview{tag,spaceing,ls}.nov_distbtwnfix; ...
                distbtwnfix(ind(1),:)];
            allview{tag,spaceing,ls}.nov_fixduration = [allview{tag,spaceing,ls}.nov_fixduration; ...
                fixduration(ind(1),:)];
            allview{tag,spaceing,ls}.nov_sacangle = [allview{tag,spaceing,ls}.nov_sacangle; ...
                sacangle(ind(1),:)];
            allview{tag,spaceing,ls}.nov_sacdist = [allview{tag,spaceing,ls}.nov_sacdist; ...
                sacdist(ind(1),:)];
            allview{tag,spaceing,ls}.nov_sacamplitude = [allview{tag,spaceing,ls}.nov_sacamplitude; ...
                sacamplitude(ind(1),:)];
            allview{tag,spaceing,ls}.nov_sacduration = [allview{tag,spaceing,ls}.nov_sacduration; ...
                sacduration(ind(1),:)];
            allview{tag,spaceing,ls}.nov_timebtwfix = [allview{tag,spaceing,ls}.nov_timebtwfix; ...
                timebtwfix(ind(1),:)];
            
            fixations = fixationstats{ind(2)}.fixations;
            for f = 2:size(fixations,2);  %ignore first since it's a fixaiton on the cross hair
                xx = round(fixations(1,f));
                yy = round(imageY-fixations(2,f));
                xx(xx < 1) = 1;
                yy(yy < 1)= 1;
                xx(xx > imageX) = imageX;
                yy(yy > imageY) = imageY;
                allview{tag,spaceing,ls}.fam_densitymap(yy,xx) = ...
                    allview{tag,spaceing,ls}.fam_densitymap(yy,xx)+1;
            end
            allview{tag,spaceing,ls}.fam_anglebtwfix = [allview{tag,spaceing,ls}.fam_anglebtwfix; ...
                anglebtwfix(ind(2),:)];
            allview{tag,spaceing,ls}.fam_sacangle_2fix = [allview{tag,spaceing,ls}.fam_sacangle_2fix; ...
                sacangle_2fix(ind(2),:)];
            allview{tag,spaceing,ls}.fam_distbtwnfix = [allview{tag,spaceing,ls}.fam_distbtwnfix; ...
                distbtwnfix(ind(2),:)];
            allview{tag,spaceing,ls}.fam_fixduration = [allview{tag,spaceing,ls}.fam_fixduration; ...
                fixduration(ind(2),:)];
            allview{tag,spaceing,ls}.fam_sacangle = [allview{tag,spaceing,ls}.fam_sacangle; ...
                sacangle(ind(2),:)];
            allview{tag,spaceing,ls}.fam_sacdist = [allview{tag,spaceing,ls}.fam_sacdist; ...
                sacdist(ind(2),:)];
            allview{tag,spaceing,ls}.fam_sacamplitude = [allview{tag,spaceing,ls}.fam_sacamplitude; ...
                sacamplitude(ind(2),:)];
            allview{tag,spaceing,ls}.fam_sacduration = [allview{tag,spaceing,ls}.fam_sacduration; ...
                sacduration(ind(2),:)];
            allview{tag,spaceing,ls}.fam_timebtwfix = [allview{tag,spaceing,ls}.fam_timebtwfix; ...
                timebtwfix(ind(2),:)];
        elseif length(ind) == 1; %if familiar image was not show fill with Nans to preserve structure
            fixations = fixationstats{ind(1)}.fixations;
            for f = 2:size(fixations,2); %ignore first since it's a fixaiton on the cross hair
                xx = round(fixations(1,f));
                yy = round(imageY-fixations(2,f));
                xx(xx < 1) = 1;
                yy(yy < 1)= 1;
                xx(xx > imageX) = imageX;
                yy(yy > imageY) = imageY;
                allview{tag,spaceing,ls}.nov_densitymap(yy,xx) = ...
                    allview{tag,spaceing,ls}.nov_densitymap(yy,xx)+1;
            end
            allview{tag,spaceing,ls}.nov_anglebtwfix = [allview{tag,spaceing,ls}.nov_anglebtwfix; ...
                anglebtwfix(ind(1),:)];
            allview{tag,spaceing,ls}.nov_sacangle_2fix = [allview{tag,spaceing,ls}.nov_sacangle_2fix; ...
                sacangle_2fix(ind(1),:)];
            allview{tag,spaceing,ls}.nov_distbtwnfix = [allview{tag,spaceing,ls}.nov_distbtwnfix; ...
                distbtwnfix(ind(1),:)];
            allview{tag,spaceing,ls}.nov_fixduration = [allview{tag,spaceing,ls}.nov_fixduration; ...
                fixduration(ind(1),:)];
            allview{tag,spaceing,ls}.nov_sacangle = [allview{tag,spaceing,ls}.nov_sacangle; ...
                sacangle(ind(1),:)];
            allview{tag,spaceing,ls}.nov_sacdist = [allview{tag,spaceing,ls}.nov_sacdist; ...
                sacdist(ind(1),:)];
            allview{tag,spaceing,ls}.nov_sacamplitude = [allview{tag,spaceing,ls}.nov_sacamplitude; ...
                sacamplitude(ind(1),:)];
            allview{tag,spaceing,ls}.nov_sacduration = [allview{tag,spaceing,ls}.nov_sacduration; ...
                sacduration(ind(1),:)];
            allview{tag,spaceing,ls}.nov_timebtwfix = [allview{tag,spaceing,ls}.nov_timebtwfix; ...
                timebtwfix(ind(1),:)];
            
            allview{tag,spaceing,ls}.fam_anglebtwfix = [allview{tag,spaceing,ls}.fam_anglebtwfix; NaN(1,75)];
            allview{tag,spaceing,ls}.fam_sacangle_2fix = [allview{tag,spaceing,ls}.fam_sacangle_2fix; NaN(1,75)];
            allview{tag,spaceing,ls}.fam_distbtwnfix = [allview{tag,spaceing,ls}.fam_distbtwnfix; NaN(1,75)];
            allview{tag,spaceing,ls}.fam_fixduration = [allview{tag,spaceing,ls}.fam_fixduration; NaN(1,75)];
            allview{tag,spaceing,ls}.fam_sacangle = [allview{tag,spaceing,ls}.fam_sacangle; NaN(1,75)];
            allview{tag,spaceing,ls}.fam_sacdist = [allview{tag,spaceing,ls}.fam_sacdist; NaN(1,75)];
            allview{tag,spaceing,ls}.fam_sacamplitude = [allview{tag,spaceing,ls}.fam_sacamplitude; NaN(1,75)];
            allview{tag,spaceing,ls}.fam_sacduration = [allview{tag,spaceing,ls}.fam_sacduration; NaN(1,75)];
            allview{tag,spaceing,ls}.fam_timebtwfix = [allview{tag,spaceing,ls}.fam_timebtwfix; NaN(1,75)];
        else %if image never show fill with NaNs to preserve structure
            allview{tag,spaceing,ls}.nov_anglebtwfix = [allview{tag,spaceing,ls}.nov_anglebtwfix; NaN(1,75)];
            allview{tag,spaceing,ls}.nov_sacangle_2fix = [allview{tag,spaceing,ls}.nov_sacangle_2fix; NaN(1,75)];
            allview{tag,spaceing,ls}.nov_distbtwnfix = [allview{tag,spaceing,ls}.nov_distbtwnfix; NaN(1,75)];
            allview{tag,spaceing,ls}.nov_fixduration = [allview{tag,spaceing,ls}.nov_fixduration; NaN(1,75)];
            allview{tag,spaceing,ls}.nov_sacangle = [allview{tag,spaceing,ls}.nov_sacangle; NaN(1,75)];
            allview{tag,spaceing,ls}.nov_sacdist = [allview{tag,spaceing,ls}.nov_sacdist; NaN(1,75)];
            allview{tag,spaceing,ls}.nov_sacamplitude = [allview{tag,spaceing,ls}.nov_sacamplitude; NaN(1,75)];
            allview{tag,spaceing,ls}.nov_sacduration = [allview{tag,spaceing,ls}.nov_sacduration; NaN(1,75)];
            allview{tag,spaceing,ls}.nov_timebtwfix = [allview{tag,spaceing,ls}.nov_timebtwfix; NaN(1,75)];
            
            allview{tag,spaceing,ls}.fam_anglebtwfix = [allview{tag,spaceing,ls}.fam_anglebtwfix; NaN(1,75)];
            allview{tag,spaceing,ls}.fam_sacangle_2fix = [allview{tag,spaceing,ls}.fam_sacangle_2fix; NaN(1,75)];
            allview{tag,spaceing,ls}.fam_distbtwnfix = [allview{tag,spaceing,ls}.fam_distbtwnfix; NaN(1,75)];
            allview{tag,spaceing,ls}.fam_fixduration = [allview{tag,spaceing,ls}.fam_fixduration; NaN(1,75)];
            allview{tag,spaceing,ls}.fam_sacangle = [allview{tag,spaceing,ls}.fam_sacangle; NaN(1,75)];
            allview{tag,spaceing,ls}.fam_sacdist = [allview{tag,spaceing,ls}.fam_sacdist; NaN(1,75)];
            allview{tag,spaceing,ls}.fam_sacamplitude = [allview{tag,spaceing,ls}.fam_sacamplitude; NaN(1,75)];
            allview{tag,spaceing,ls}.fam_sacduration = [allview{tag,spaceing,ls}.fam_sacduration; NaN(1,75)];
            allview{tag,spaceing,ls}.fam_timebtwfix = [allview{tag,spaceing,ls}.fam_timebtwfix; NaN(1,75)];
        end
    end
end

% Combine all novel data across monkeys spaceing and lesion status
nov_densitymap = zeros(imageY,imageX); %PDF of fixation locations
nov_anglebtwfix = []; %angle between fixation locations
nov_sacangle_2fix = []; %saccade angle entering a fixation
nov_distbtwnfix = [];%distance between fixations
nov_fixduration = [];% fixaiton durations
nov_sacangle = [];%saccade angle leaving a fixation
nov_sacdist = [];%saccade arc length
nov_sacamplitude = [];%saccade amplitude
nov_sacduration = [];%saccade duration
nov_timebtwfix = [];%time between fixations
for t = 1:length(tags)
    for space = 1:5
        for ls = 1:2
            nov_densitymap = nov_densitymap + allview{t,space,ls}.nov_densitymap;
            nov_anglebtwfix = [nov_anglebtwfix; allview{t,space,ls}.nov_anglebtwfix];
            nov_sacangle_2fix = [nov_sacangle_2fix; allview{t,space,ls}.nov_sacangle_2fix];
            nov_distbtwnfix = [nov_distbtwnfix; allview{t,space,ls}.nov_distbtwnfix];
            nov_fixduration = [nov_fixduration; allview{t,space,ls}.nov_fixduration];
            nov_sacangle = [nov_sacangle; allview{t,space,ls}.nov_sacangle];
            nov_sacdist = [nov_sacdist; allview{t,space,ls}.nov_sacdist];
            nov_sacamplitude = [nov_sacamplitude; allview{t,space,ls}.nov_sacamplitude];
            nov_sacduration = [nov_sacduration; allview{t,space,ls}.nov_sacduration];
            nov_timebtwfix = [nov_timebtwfix; allview{t,space,ls}.nov_timebtwfix];
        end
    end
end

f = fspecial('gaussian',[256,256],24); %1 dva 2D gaussian filter
nov_densitymap = imfilter(nov_densitymap,f);
figure
imagesc(nov_densitymap)
title('PDF of fixation locations during novel presentations-all monkeys and spacings')
axis off

nov_anglebtwfix = nov_anglebtwfix(1:end);
nov_anglebtwfix(isnan(nov_anglebtwfix)) = [];
[probanglebtwfix] = hist(nov_anglebtwfix,360);
probanglebtwfix = [probanglebtwfix(36:-1:1) probanglebtwfix probanglebtwfix(end:-1:end-36)];
probanglebtwfix = filtfilt(1/18*ones(1,18),1,probanglebtwfix);
probanglebtwfix = probanglebtwfix(37:end-37);
probanglebtwfix = probanglebtwfix/sum(probanglebtwfix);
probanglebtwfix = [probanglebtwfix probanglebtwfix(1)];
polar(n,100*probanglebtwfix); %n is 0 to 2pi
title('PDF of angle between fixations')

nov_sacangle = nov_sacangle(1:end);
nov_sacangle(isnan(nov_sacangle)) = [];
[probsacangle] = hist(nov_sacangle,360);
probsacangle = [probsacangle(36:-1:1) probsacangle probsacangle(end:-1:end-36)];
probsacangle = filtfilt(1/6*ones(1,6),1,probsacangle);
probsacangle = probsacangle(37:end-37);
probsacangle = probsacangle/sum(probsacangle);
probsacangle = [probsacangle probsacangle(1)];
polar(n,100*probsacangle); %n is 0 to 2pi
title('PDF of saccade angles leaving a fixation')

nov_sacangle_2fix = nov_sacangle_2fix(1:end);
nov_sacangle_2fix(isnan(nov_sacangle_2fix)) = [];
[probsacangle_2fix] = hist(nov_sacangle_2fix,360);
probsacangle_2fix = [probsacangle_2fix(36:-1:1) probsacangle_2fix probsacangle_2fix(end:-1:end-36)];
probsacangle_2fix = filtfilt(1/6*ones(1,6),1,probsacangle_2fix);
probsacangle_2fix = probsacangle_2fix(37:end-37);
probsacangle_2fix = probsacangle_2fix/sum(probsacangle_2fix);
probsacangle_2fix = [probsacangle_2fix probsacangle_2fix(1)];
polar(n,100*probsacangle_2fix); %n is 0 to 2pi
title('PDF of saccade angles entering fixations')

mean_novdistbtwnfix = nanmean(nov_distbtwnfix);
figure
plot(nanmean(nov_distbtwnfix(:,1:25))/24)
box off
ylabel('Distance (dva)')
xlabel('Fixation Number')
title('Distance between fixations')

novdistfix = nov_distbtwnfix(:,1:25)/24;
figure
hist(novdistfix(1:end),25)
xlabel('Distance (dva)')
ylabel('count')
title('PDF of distance between fixations')

figure
plot(samprate*nanmean(nov_fixduration(:,1:35)))
xlabel('Fixation number')
ylabel('Duration (ms)')
box off

novfixdur = samprate*nov_fixduration(:,1:25);
figure
hist(novfixdur(1:end),50)
xlabel('Fixation duration (ms)')
ylabel('count')
box off

figure
plot(nanmean(nov_sacdist(:,1:35))/24)
xlabel('Saccade Number')
ylabel('Saccade Arc Length (dva)')
box off

figure
plot(nanmean(nov_sacamplitude(:,1:35)/24))
xlabel('Saccade Number')
ylabel('Saccade Amplitude (dva)')
box off

novsacamp = nov_sacamplitude(:,1:25)/24;
figure
hist(novsacamp(1:end),50)
xlabel('Saccade Amplitude (dva)')
ylabel('Count')

novsacdur = samprate*nov_sacduration(:,1:25);
novsacdur(novsacdur > 150) = [];
figure
hist(novsacdur(1:end),25)
xlabel('Saccade Duration (ms)')
ylabel('Count')
box off

novsacrate = 1000./(samprate*nov_timebtwfix(:,1:25));
figure
hist(novsacrate(1:end),50)
xlabel('Instantatneious saccade rate (Hz)')
ylabel('Count')
box off

fixationduration = cell(2,5,2); %row is novel vs familiar and col is healthy vs lesion
saccadeamplitudes = cell(2,5,2);%row is novel vs familiar and col is healthy vs lesion

for t = 1:length(tags);
    for space = 1:5;
        for ls = 1:2
            fixationduration{1,space,ls} = [fixationduration{1,space,ls}; ...
                allview{t,space,ls}.nov_fixduration];
            fixationduration{2,space,ls} = [fixationduration{1,space,ls}; ...
                allview{t,space,ls}.fam_fixduration];
            
            saccadeamplitudes{1,space,ls} = [saccadeamplitudes{1,space,ls}; ...
                allview{t,space,ls}.nov_fixduration];
            saccadeamplitudes{2,space,ls} = [saccadeamplitudes{1,space,ls}; ...
                allview{t,space,ls}.fam_fixduration];
        end
    end
end

figure
hold all
for novelty = 1:2;
    for space = 1:5;
        plot(samprate*nanmedian(fixationduration{novelty,space,1}(:,1:25)));
    end
end
legend('Novel 0','Novel 2','Novel 4','Novel 8','Novel 16','Familiar 0',...
    'Familiar 2','Familiar 4','Familiar 8','Familiar 16')
xlabel('Fixation Number')
ylabel('Fixation Duration (ms)')
title('Healthy Monkeys: Fixation duration by fixation number novel vs familiar presentations')


figure
for space = 1:5
    subplot(3,2,space)
    hold on
    novfixdur = samprate*fixationduration{1,space,1};
     numpoints = sum(~isnan(novfixdur(1:end)));
    errorbar(1,nanmean(novfixdur(1:end)),nanstd(novfixdur(1:end))/sqrt(numpoints),'b')
    famfixdur = samprate*fixationduration{2,space,1};
    numpoints = sum(~isnan(famfixdur(1:end)));
    errorbar(2,nanmean(famfixdur(1:end)),nanstd(famfixdur(1:end))/sqrt(numpoints),'r')
    hold off
    title(['Spaceing of ' num2str(2^(space-1))])
end
title('Healthy Monkeys: Fixation duration across all fixations novel vs familiar presentations')

figure
hold all
for novelty = 1:2;
    for space = 1:5;
        plot(samprate*nanmedian(fixationduration{novelty,space,2}(:,1:25)));
    end
end
legend('Novel 0','Novel 2','Novel 4','Novel 8','Novel 16','Familiar 0',...
    'Familiar 2','Familiar 4','Familiar 8','Familiar 16')
xlabel('Fixation Number')
ylabel('Fixation Duration (ms)')
title('Lesion Monkeys: Fixation duration by fixation number novel vs familiar presentations')


figure
for space = 1:5
    subplot(3,2,space)
    hold on
    novfixdur = samprate*fixationduration{1,space,2};
     numpoints = sum(~isnan(novfixdur(1:end)));
    errorbar(1,nanmean(novfixdur(1:end)),nanstd(novfixdur(1:end))/sqrt(numpoints),'b')
    famfixdur = samprate*fixationduration{2,space,2};
    numpoints = sum(~isnan(famfixdur(1:end)));
    errorbar(2,nanmean(famfixdur(1:end)),nanstd(famfixdur(1:end))/sqrt(numpoints),'r')
    hold off
    title(['Spaceing of ' num2str(2^(space-1))])
end
subtitle('Lesion Monkeys: Fixation duration across all fixations novel vs familiar presentations')