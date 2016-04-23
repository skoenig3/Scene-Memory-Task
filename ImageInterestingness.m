SMT_dir = 'C:\Users\seth.koenig\Documents\MATLAB\Scene Memory Task\'; %parent directory
cd(SMT_dir)

imageX = 800; %horizontal pixel size of images
imageY = 600; %vertical pixel size of images

% sobel filters detect edges!
sobelx = [1     2   1;
    0     0   0;
    -1    -2  -1;];

sobely = [1     2   1;
    0     0   0;
    -1    -2  -1;];

load('allcortexfiles.mat'); %easily accessible data structure from above

entropyvalues = NaN(2,625);
edgevalues = NaN(2,625);
lookingtime = NaN(2,625);
count = [1,1];
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
    
    
    numrpt = size(event_arr,2);
    valrptcnt = 0;
    clear per clrchgind
    for rptlop = 1:numrpt
        if itmlist(event_arr((find(event_arr(:,rptlop)>1000,1,'last')),rptlop)-1000) <= 79
            if size(find(event_arr(:,rptlop) == 200)) ~=0
                perbegind = find(event_arr(:,rptlop) == 24);%was originally 23, changed this and begtimdum line below to optimize
                perendind = find(event_arr(:,rptlop) == 24);
                if length( perbegind) > 1
                    perbegind = perbegind(2);
                    perendind = perendind(2);
                end
                cndnumind = find(event_arr(:,rptlop) >= 1000 & event_arr(:,rptlop) <=2000);
                blknumind = find(event_arr(:,rptlop) >=500 & event_arr(:,rptlop) <=999);
                begtimdum = time_arr(perbegind,rptlop)-100;
                endtimdum = time_arr(perendind,rptlop);
                if endtimdum > begtimdum
                    valrptcnt = valrptcnt + 1;
                    clrchgind(valrptcnt)=rptlop;
                    per(valrptcnt).begsmpind = begtimdum;
                    per(valrptcnt).endsmpind = endtimdum;
                    per(valrptcnt).begpos = 1;
                    per(valrptcnt).cnd = event_arr(cndnumind,rptlop);
                    per(valrptcnt).blk = event_arr(blknumind,rptlop);
                    per(valrptcnt).allval = event_arr(:,rptlop);
                    per(valrptcnt).alltim = time_arr(:,rptlop);
                    per(valrptcnt).event = rptlop;
                end
            end
        end
    end
    
    samprate=5;
    
    per = per(3:end);
    %Don't keep first 2 calibration pionts these are for offset correction at
    %start of task
    clear cnd
    numrpt = size(per,2);
    cnd = zeros(1,numrpt);
    for rptlop = 1:numrpt
        cnd(rptlop)=per(rptlop).cnd;
    end
    
    %---For Calibration with Eye tracking data with cp2tform---%
    % Create structures x and y of the corresponding average eye data for each trial
    % instance (l) of each condition (k)
    spacex = [-12,-6,0,6,12];
    spacey = [-8,-4,0,4,8];
    x = cell(5,5);
    y = cell(5,5);
    control = NaN(length(cnd),2);
    % clr = ['rgbmkrgbmkrgbmkrgbmkrgbmk'];
    % figure
    % hold on
    for k = 1:length(cnd)
        C = textscan(itmfil(itmlist(cnd(k)-1000)+6,:),'%d');
        control(k,:) = C{1}(4:5)';
        
        xi = find(C{1}(4) == spacex);
        yi = find(C{1}(5) == spacey);
        eyeind = floor(((per(k).begsmpind-1000)/samprate)*2):(floor((per(k).endsmpind-1000)/samprate))*2;
        evenind = eyeind(logical(~rem(eyeind,2)));
        oddind =  eyeind(logical(rem(eyeind,2)));
        x{xi,yi} = [x{xi,yi} mean(eog_arr(oddind,per(k).event))];
        y{xi,yi} = [y{xi,yi} mean(eog_arr(evenind,per(k).event))];
        if strcmpi(cortexfile(end-9:end),'PW131028.2') %calibration issue with this file so fixing eye traces here
            if cnd(k)-1000 < 125
                y{xi,yi} = [y{xi,yi} mean(eog_arr(evenind,per(k).event))+200];
            end
        elseif strcmpi(cortexfile(end-9:end),'PW131101.2') %calibration issue with this file so fixing eye traces here
            if cnd(k)-1000 >= 127
                y{xi,yi} = [y{xi,yi} mean(eog_arr(evenind,per(k).event))-200];
            end
        end
        
        %     plot(mean(eog_arr(oddind,per(k).event)),mean(eog_arr(evenind,per(k).event)),[clr(xi*yi) '+'])
    end
    % title(['Calibration transformation for ' cortexfile(end-9:end)])
    
    %Test for errors%
    count = zeros(5,5);
    for xi = 1:length(spacex);
        for yi = 1:length(spacey);
            count(yi,xi) = sum(control(:,1) == spacex(xi) & control(:,2) == spacey(yi));
        end
    end
    if any(count ~= 10);
        disp('Calibration trial analysis incomplete or error')
        disp('Check number of calibration pionts or task not finished')
    end
    
    clear meanx meany
    for k=1:numel(x)
        xss = x{k};
        low = mean(xss)-std(xss);
        high = mean(xss)+std(xss);
        xss(xss < low) = [];
        xss(xss > high) = [];
        meanx(k)=median(xss);
    end
    for k=1:numel(y)
        yss = y{k};
        low = mean(yss)-std(yss);
        high = mean(yss)+std(yss);
        yss(yss < low) = [];
        yss(yss > high) = [];
        meany(k)=median(y{k});
    end
    
    controlx = [];
    controly = [];
    for i = 1:length(spacex);
        for ii = 1:length(spacey);
            controly = [controly spacey(i)];
            controlx = [controlx spacex(ii)];
        end
    end
    
    tform = cp2tform([controlx' controly'], [meanx' meany'],'polynomial',4);
    tform.forward_fcn = tform.inverse_fcn;
    
    figure
    hold on
    for i = 1:length(controlx);
        plot(controlx(i),controly(i),'r+')
        [x,y] = tformfwd(tform,meanx(i),meany(i));
        plot(x,y,'*b')
    end
    title(['Calibration transformation for ' cortexfile(end-9:end)])
    
    numrpt = size(event_arr,2);
    valrptcnt = 0;
    clear per vpcind
    new_eog_arr=[];
    for rptlop = 1:numrpt
        itmlist(event_arr((find(event_arr(:,rptlop)>1000,1,'last')),rptlop)-1000);
        if itmlist(event_arr((find(event_arr(:,rptlop)>1000,1,'last')),rptlop)-1000) >= 80
            if size(find(event_arr(:,rptlop) == 200)) ~=0
                perbegind = find(event_arr(:,rptlop) == 23,1,'first');
                perendind = find(event_arr(:,rptlop) == 24,1,'first');
                cndnumind = find(event_arr(:,rptlop) >= 1000 & event_arr(:,rptlop) <=2000);
                blknumind = find(event_arr(:,rptlop) >=500 & event_arr(:,rptlop) <=999);
                begtimdum = time_arr(perbegind,rptlop);
                endtimdum = time_arr(perendind,rptlop);
                if endtimdum > begtimdum
                    valrptcnt = valrptcnt + 1;
                    vpcind(valrptcnt)=rptlop;
                    per(valrptcnt).begsmpind = begtimdum;
                    per(valrptcnt).endsmpind = endtimdum;
                    per(valrptcnt).begpos = 1;
                    per(valrptcnt).cnd = event_arr(cndnumind,rptlop);
                    per(valrptcnt).blk = event_arr(blknumind,rptlop);
                    per(valrptcnt).allval = event_arr(:,rptlop);
                    per(valrptcnt).alltim = time_arr(:,rptlop);
                    new_eog_arr=cat(2,new_eog_arr,eog_arr(:,rptlop));
                end
            end
        end
    end
    
    %---get eye data for only when fixation cross or picture is displayed---%
    eyedat = cell(1,length(per));
    for trlop=1:size(per,2)
        trleog=new_eog_arr(~isnan(new_eog_arr(:,trlop)),trlop); % eog for this trial
        horeog=trleog(1:2:size(trleog,1)); % horizontal eye dat
        vrteog=trleog(2:2:size(trleog,1)); % vertical eye dat
        picstart=per(trlop).alltim(find(per(trlop).allval==23,1,'first'))-per(trlop).alltim(find(per(trlop).allval==100)); % picture start time relative to eye scan start
        picend=per(trlop).alltim(find(per(trlop).allval==24,1,'first'))-per(trlop).alltim(find(per(trlop).allval==100)); % picture end time relative to eye scan start
        
        if picend/5<=length(horeog)
            eyedat{trlop}(1,:) = (horeog(round(picstart/5):floor(picend/5)));
            eyedat{trlop}(2,:) = (vrteog(round(picstart/5):floor(picend/5)));
            if strcmpi(cortexfile(end-9:end),'PW131028.2') %calibration issue with this file so fixing eye traces here
                if per(trlop).cnd-1000 < 125
                    eyedat{trlop}(2,:) = (vrteog(round(picstart/5):floor(picend/5)))+200;
                end
            elseif strcmpi(cortexfile(end-9:end),'PW131101.2') %calibration issue with this file so fixing eye traces here
                if per(trlop).cnd-1000 >= 127
                    eyedat{trlop}(2,:) = (vrteog(round(picstart/5):floor(picend/5)))-200;
                end
            end
        else
            eyedat{trlop}(1,:)=nan;
            eyedat{trlop}(2,:)=nan;
        end
    end
    
    cnd=[];
    numrpt = size(per,2);
    for rptlop = 1:numrpt
        cnd(rptlop)=per(rptlop).cnd;
    end
    
    %---Recalibrate and automatically scale eye data---%
    for eye = 1:length(eyedat)
        x = eyedat{eye}(1,:);
        y = eyedat{eye}(2,:);
        [x,y] = tformfwd(tform,x,y);
        eyedat{eye} = [x;y];
    end
    
    images = zeros(length(eyedat),1);
    for i = 1:length(cnd);
        C = textscan(itmfil(itmlist(cnd(i)-1000)+6,:),'%s');
        img = C{1}{end};
        slash = strfind(img,'\');
        img = img(slash(2)+1:end-4);
        images(i) = str2double(img);
    end
    
    for im = 1:25;
        ind = find(img == im);
        if ~isempty(ind);
            
            imge = imread([image_dir num2str(im) '.bmp']);
            imge = rgb2gray(imge);
            ent = entropy(imge);%pixel intesnity entropy
            xedges = imfilter(imge,sobelx);
            yedges = imfilter(imge,sobely);
           	edge = mean2(xedges+yedges); %edgineess
            
            if length(ind) == 1;
                entropyvalues(1,count(1)) =ent;
                edgevalues(1,count(1));
                lookingtime(1,count(1)) = 5*size(eyedat{ind(2)},2);
                count(1) = count(1)+1;
            end
            if length(ind) == 2
                
                entropyvalues(2,count(2));
                edgevalues(2,count(2));
                lookingtime(2,count(2)) = 5*size(eyedat{ind(2)},2);
                count(2) = count(2)+1;
            end
        end
    end
    
    
    for i = 1:size(a,1);
        index = strfind(a(i,:),'jpg');
        if isempty(index);
            index = strfind(a(i,:),'jpeg');
        end
        if ~isempty(index);
            img = imread(a(i,:));
            img = rgb2gray(img);
            entropyvalues(i) = entropy(img);%pixel intesnity entropy
            xedges = imfilter(img,sobelx);
            yedges = imfilter(img,sobely);
            edgevalues(i) = mean2(xedges+yedges); %edgineess
        end
        
    end
    
end
