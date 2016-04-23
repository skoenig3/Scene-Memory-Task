%% Generate images sets for Scene Memory Task (SMT)
% [1] Grab Potential SMT Pictures from Flickr
% [2] Remove images that are too small or are grayscale
% [3] remove images low in image entropy and edginess (below 50% on each)
% [4] Rename files if this helps
% [5] Resize images
% [6] Get Salience maps & Salience Map entropy if this helps
% [7] Rename files in set folder and change to bmp if not already bmp
%% [1] Grab Potential SMT Pictures from Flickr
% function [] = getflickr(numpics)
% Mike Jutras
% Nathan Killian
% last updated 6/19/2012 nk, fixed url names and now grabs hi-res (_b) photos

% may be a max of 500 per day, 3500 for 7 days
% if photo is already in the folder it will be overwritten

% number of pictures desired
numpics=3500;

saveDir = 'C:\Users\seth.koenig\Documents\MATLAB\Scene Memory Task\Flickr Images\Flickr Images';


l=0;%number of pics grabbed
for m=1:1000% number of reloads of the page
    s=urlread('http://www.flickr.com/explore/interesting/7days/');% (re)load the site source code
    
    clear siteind
    hind=strfind(s,'h');
    i=1;
    for k=1:length(hind)
        if ((hind(k)+6)<length(s))
            if strcmp(s(hind(k):hind(k)+6),'http://')
                firstquote=strfind(s(hind(k):end),'"');
                siteind{i}=s(hind(k):(firstquote(1)+hind(k)-2));
                i=i+1;
            end
        end
    end
    
    clear imageind
    i=1;
    for k=1:size(siteind,2)
%         siteind{k}
        if size(siteind{k},2)>30 %fixme: make generic 1-9
            if strcmp(siteind{k}(1:30),'http://farm6.staticflickr.com/')
                imageind{i}=siteind{k}(1:63); %63 char limit added 30-Jul-2009 due to error in flickr html source
                i=i+1;
            elseif strcmp(siteind{k}(1:30),'http://farm8.staticflickr.com/')
                imageind{i}=siteind{k}(1:63); %additional image server added 1/15/10
                i=i+1;
            elseif strcmp(siteind{k}(1:30),'http://farm9.staticflickr.com/')
                imageind{i}=siteind{k}(1:63); %additional image server added 1/15/10
                i=i+1;
            elseif strcmp(siteind{k}(1:30),'http://farm5.staticflickr.com/')
                imageind{i}=siteind{k}(1:63); %additional image server added 1/15/10
                i=i+1;
            elseif strcmp(siteind{k}(1:30),'http://farm4.staticflickr.com/')
                imageind{i}=siteind{k}(1:63); %additional image server added 1/15/10
                i=i+1;
            elseif strcmp(siteind{k}(1:30),'http://farm3.staticflickr.com/')
                imageind{i}=siteind{k}(1:63); %additional image server added 1/15/10
                i=i+1;
            end
        end
    end
    
    clear x
    for k=1:length(imageind)
        imageind{k} = [imageind{k}(1:end-5) 'b.jpg'];% b suffix is higher resolution, h appears to be higher res, but it looks like not all photos have an h type
        x{k}=imread(imageind{k});
        %     figure;image(x{k})
    end
    
    for k=1:length(x)
        imtitle=[saveDir imageind{k}(find(double(imageind{k})==47,1,'last')+1:end)];
        if exist(imtitle,'file')~=2
            %             imwrite(x{k},imtitle,'jpg','Bitdepth',12,'Mode','lossless','Quality',100);
            imwrite(x{k},[imtitle(1:end-4) '.bmp'],'BMP');
            l=l+1;
        end
        if l==numpics
            break
        end
    end
    
    if l==numpics
        break
    end
    
end
%% [2] Remove images that are too small or are grayscale
imgdir = 'C:\Users\seth.koenig\Documents\MATLAB\Scene Memory Task\Flickr Images';
cd(imgdir)
smalldir = 'Images that are too small';
mkdir(smalldir);
graydir = 'Gray scale images';
mkdir(graydir);
a = ls;

toosmallimages = [];
grayimages = [];
for i = 1:size(a,1);
    index = strfind(a(i,:),'bmp');
    if isempty(index);
        index = strfind(a(i,:),'jpg');
        if isempty(index)
             index = strfind(a(i,:),'jpeg');
        end
    end
    if ~isempty(index);
        img = imread(a(i,:));
        if size(img,3) == 1
            grayimages= [grayimages i];
        elseif all(all(img(:,:,1) == img(:,:,2))) ...
                || all(all(img(:,:,2) == img(:,:,3))) || all(all(img(:,:,1) == img(:,:,3)))
            grayimages= [grayimages i];
        elseif size(img,1) < 600 || size(img,2) < 800
            toosmallimages = [toosmallimages i];
        end
    end
end

for ti = 1:length(toosmallimages);
    movefile(a(toosmallimages(ti),:),smalldir)
end
for gi = 1:length(grayimages);
    movefile(a(grayimages(gi),:),graydir);
end

%% [3] remove images low in image entropy and edginess (below 50% on each)
clear, clc
imgdir = 'C:\Users\seth.koenig\Documents\MATLAB\Scene Memory Task\Flickr Images';
cd(imgdir)
a = ls;

% sobel filters detect edges!
sobelx = [1     2   1;
    0     0   0;
    -1    -2  -1;];

sobely = [1     2   1;
    0     0   0;
    -1    -2  -1;];

entropyvalues = NaN(1,size(a,1));
edgevalues = NaN(1,size(a,1));
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

entropythresh = nanmean(entropyvalues);
edgethresh = nanmean(edgevalues);
gonogo = zeros(1,size(a,1));
for i = 1:size(a,1);
    if entropyvalues(i) > entropythresh && edgevalues(i) > edgethresh
        gonogo(i) = 1;
    end
end

gooddir = 'Good Images';
mkdir(gooddir);


for i = 1:size(a,1);
    if gonogo(i) == 1
        movefile(a(i,:),gooddir);
    end
end
%% [4] Rename files if this helps
% can skip section not essential
clear, clc
imgdir = 'C:\Users\seth.koenig\Documents\MATLAB\Scene Memory Task\Flickr Images\Good Images';
cd(imgdir)
a = ls;

for aa = 3:size(a,1);
    movefile(a(aa,:),[num2str(aa-2) '.bmp'])
end
%% [5] Resize images
% assumes no .mat files in the folder
% crops image borders so takes central 800 x 600 pixels
clear, clc
imgdir = 'C:\Users\seth.koenig\Documents\MATLAB\Scene Memory Task\Flickr Images\Good Images';
cd(imgdir)
a = ls;
for aa = 3:size(a,1);
    img = imread(a(aa,:));
    if all(size(img) == [600,800,3]);
        imwrite(img,a(aa,:),'bmp');
    else
        if size(img,1) ~= 600
            extra = size(img,1)-600;
            if rem(extra,2) ~= 0;
                top = ceil(extra/2);
                bottom = floor(extra/2);
            else
                top = extra/2;
                bottom = extra/2;
            end
            img = img(top:end-bottom-1,:,:);
        end
        if size(img,2) ~= 800
            extra = size(img,2)-800;
            if rem(extra,2) ~= 0;
                left = ceil(extra/2);
                right = floor(extra/2);
            else
                left = extra/2;
                right = extra/2;
            end
            img = img(:,left:end-right-1,:);
        end
        imwrite(img,a(aa,:),'bmp');
    end
end
%% [6] Get Salience maps & Salience Map entropy if this helps
% Very high entropy in salience maps and well pretty much normally distributed with small std
% Can ignore section
clear, clc
imgdir = 'C:\Users\seth.koenig\Documents\MATLAB\Scene Memory Task\Flickr Images\Good Images';
cd(imgdir)
a = ls;
for aa = 3:size(a,1);
    getSalienceMap(a(aa,:));
end

clear, clc
imgdir = 'C:\Users\seth.koenig\Documents\MATLAB\Scene Memory Task\Flickr Images\Good Images';
cd(imgdir)
a = what;
values = NaN(1,500);
for i = 1:length(a.mat);
    load(a.mat{i},'fullmap');
    values(i) = entropy(fullmap);
end
values(isnan(values)) = [];

% Display salience maps to see which might be good
clear, clc
imgdir = 'C:\Users\seth.koenig\Documents\MATLAB\Scene Memory Task\Flickr Images\Good Images';
cd(imgdir)
a = ls;
for aa = 50:size(a,1);
    load([num2str(aa-2) '-saliencemap.mat'],'fullmap');
    figure,imagesc(fullmap);
    title(['Image # ' num2str(aa-2)])
    axis off
end

%% [7] Rename files in set folder and change to bmp if not already bmp
curdir = 'C:\Users\seth.koenig\Desktop\';
newdir = 'C:\Users\seth.koenig\Desktop\SMT005\';
a = ls;

cd(curdir)
mkdir(newdir)
for i = 3:size(a,1);
    img = imread(a(i,:));
    imwrite(img,[newdir num2str(i-2) '.bmp'],'bmp');
end