function getSMTeyedat(cortexfile,ITMFile,CNDFile,imageX,imageY)
% extracts raw Shift task (SFT) eye data from a cortex files, calibrates
% the eye data, and then uses Cluster Fix to extract fixations and saccades

itmfil=[];
h =fopen(ITMFile);
tline = 1;
while tline > 0
    tline = fgetl(h);
    if ~isempty(itmfil)
        if length(tline)>size(itmfil,2)
            tline=tline(1:size(itmfil,2));
        end
    end
    tline = [tline ones(1,(size(itmfil,2)-length(tline)))*char(32)];
    if ischar(tline)
        itmfil=[itmfil; tline];
    else
        break
    end
end
fclose(h);

cndfil=[];
h=fopen(CNDFile);
tline = 1;
while tline > 0
    tline = fgetl(h);
    if ~isempty(cndfil)
        if length(tline)>size(cndfil,2)
            tline=tline(1:size(cndfil,2));
        end
    end
    tline = [tline ones(1,(size(cndfil,2)-length(tline)))*char(32)];
    if ischar(tline)
        cndfil=[cndfil; tline];
    else
        break
    end
end
fclose(h);

itmlist = zeros(size(cndfil,1)-1,1);
for i = 2:size(cndfil,1);
    str = textscan(cndfil(i,:),'%d');
    itmlist(i-1) = str{1}(end);
end

[time_arr,event_arr,eog_arr,epp_arr,header,trialcount]  = get_ALLdata(cortexfile);

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


for i = 1:size(eyedat,2);
    x = 24*eyedat{i}(1,:)+imageX/2;
    y = 24*eyedat{i}(2,:)+imageY/2;
    %     if length(x) > 1520 %500 ms cross hair + 7000 ms for image + 100 ms buffer
    %        x = x(1:1520);
    %        y = y(1:1520);
    %     end
    badx = find(x < 50 | x > imageX+50); %~1 dva leave margin of error
    x(badx) = []; y(badx) = [];
    bady = find(y < -50 | y > imageY+50); %~1 dva margin of error
    x(bady) = []; y(bady) = [];
    eyedat{i} = [x;y];
    %     C = textscan(itmfil(itmlist(cnd(i)-1000)+6,:),'%s');
    %     img = C{1}{end};
    %     slash = strfind(img,'\');
    %     img = img(slash(2)+1:end-4);
    %     img = imread([pwd '\Image Sets\' ITMFile(end-9:end-4) '\' num2str(img) '.bmp']);
    %     screen_size = get(0, 'ScreenSize');
    %     figure
    %     axis off
    %     imshow(img)
    %     hold on
    %     plot(x,imageY-y,'b')
    %     set(gcf, 'Position', [0 0 screen_size(3) screen_size(4) ] );
end


images = zeros(length(eyedat),1);
for i = 1:length(cnd);
    C = textscan(itmfil(itmlist(cnd(i)-1000)+6,:),'%s');
    img = C{1}{end};
    slash = strfind(img,'\');
    img = img(slash(2)+1:end-4);
    images(i) = str2double(img);
end


[fixationstats] = ClusterFixation_Final(eyedat);
save([cortexfile(1:56) 'datafiles\' cortexfile(end-9:end) '-fixation.mat'],...
    'images','fixationstats')
end