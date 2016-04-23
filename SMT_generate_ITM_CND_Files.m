%% Create Item (.ITM) and Condition (.CND) files for Scene Memory Task (SMT)
% [1] Generate random sequences of images via brute force.
% [2] Generate Item/CND files 
% [3] Generate Item files for wide screen setup(s)

%% [1] Generate random sequences of images via brute force.
% Saves sequences at end so check so you don't overwrite!!!
% Creates sequences with spacing appart indicated in variable named
% "spaceing" with the number of novel images indicated in variable
% "numimages". There are an equal number of repeated presentations as there
% are novel presentations. 
clc,clear

spaceing = [0 2 4 8 16];
numimages = [5 5 5 5 5];
totalimages = sum(numimages);

sequences = NaN(20,totalimages*2);

success = 0;
supercount = 0;
while success < 20 && supercount < 100
    count = 0;
    nans = 1;
    supercount = supercount + 1;
    while count <= 100000 && nans ~= 0
        count = count + 1;
        images = NaN(1,totalimages*2);
        imgnum = 1;
        for s = length(spaceing):-1:1;
            maxinit = length(images)-spaceing(s)-1;
            usedind = find(~isnan(images));
            posind = 1:maxinit;
            [~,ia,~]=intersect(posind,usedind);
            posind(ia) = [];
            [~,ia,~]=intersect(posind+spaceing(s)+1,usedind);
            posind(ia) = [];
            starts = posind(randperm(length(posind)));
            if length(starts) >=numimages(s);
                starts = starts(1:numimages(s));
                for i = 1:length(starts);
                    images(starts(i)) = imgnum;
                    images(starts(i)+spaceing(s)+1) = imgnum;
                    imgnum = imgnum+1;
                end
            else
                break
            end
        end
        nans = sum(isnan(images));
        if nans == 0
            
            % check following errors
            imgcount = zeros(1,totalimages);
            for i = 1:totalimages;
                imgcount(i) = sum(images == i);
            end
            if all(imgcount == 2);
                disp('Correct number of images')
                space = zeros(1,totalimages);
                for i = 1:totalimages;
                    ind = find(images == i);
                    space(i) = diff(ind)-1;
                end
                spaces = space;
                s = length(spaceing)+1;
                while s >= 2
                    s = s -1;
                    if length(find(space == spaceing(s))) ~= numimages(s)
                        disp('error in spacing')
                        s = NaN;
                    end
                end
                if s == 1;
                    disp('spacing is good')
                    firstnan = find(isnan(sequences(:,1)));
                    firstnan = firstnan(1);
                    disp(['Sequence # ' num2str(firstnan) ' determined'])
                    sequences(firstnan,:) = images;
                    success = success + 1;
                end
            else
                disp('missing or extra images')
            end
        end
    end
end
saveas('imagesequences.mat','sequences');
%% [2] Generate Item/CND files 
% consistently but differeintly with random seed NOT for wide screen version
numimages = 25;
for ssss = 2:5
    setnum = ssss;
    rand('seed',setnum) % does not apply to set001.itm
    if setnum < 10
        set = ['SMT00' num2str(setnum)];
    else
        set = ['SMT0' num2str(setnum)];
    end
    
    fid = fopen([set '.itm'],'w+');
    
    spacingx = {'-12.0','-6.00',' 0.00',' 6.00',' 12.0'};
    spacingy = {'-8.00','-4.00',' 0.00',' 4.00',' 8.00'};
    ngrid = length(spacingx);
    
    line1 ='ITEM TYPE FILLED CENTERX CENTERY BITPAN HEIGHT WIDTH ANGLE OUTER INNER -R- -G- -B- C A ------FILENAME------';
    line2 ='  -4    1      1    0.00    0.00      0   1.00  1.00  0.00               0   0   0 x';
    line3 ='  -3   14      1    0.00    0.00      0   0.50  0.50  0.00             255 255 255 x';
    line4 ='  -2    1      1    0.00    0.00      0   0.20  0.20  0.00  0.50       100 100 100 x';
    line5 ='  -1    1      0    0.00    0.00      0   1.00  1.00  0.00             255 255 255 x';
    line6 ='   0    1      1    0.00    0.00      0   0.15  0.15  0.00              37  63  49 x';
    line7 ='   1    1      1    0.00    0.00      0   0.30  0.30  0.00  0.50       150 150 150 x';
    
    
    
    for line = 1:7
        fprintf(fid,[eval(['line' num2str(line)]) '\r\n']);
    end
    for i = 1:length(spacingx)
        for ii = 1:length(spacingy)
            num1 = num2str(1+3*ngrid*(i-1)+2*(ii-1)+ii);
            num3 = num2str(1+3*ngrid*(i-1)+2*(ii-1)+ii+1);
            num2 = num2str(1+3*ngrid*(i-1)+2*(ii-1)+ii+2);
            if str2double(num1) < 10
                header1 = ['   ' num1];
            else
                header1 = ['  ' num1];
            end
            if str2double(num2) < 10
                header2 = ['   ' num2];
            else
                header2 = ['  ' num2];
            end
            if str2double(num3) < 10
                header3 = ['   ' num3];
            else
                header3 = ['  ' num3];
            end
            str1 = [header1 '    1      1   ' spacingx{i} '   ' spacingy{ii} '      0   0.30  0.30  0.00  0.50       150 150 150 x' '\r\n'];
            str3 = [header3 '    1      1   ' spacingx{i} '   ' spacingy{ii} '      0   0.30  0.30  0.00  0.50       150 150 150 x' '\r\n'];
            str2 = [header2 '    1      1   ' spacingx{i} '   ' spacingy{ii} '      0   0.30  0.30  0.00  0.50       175 175 130 x' '\r\n'];
            fprintf(fid,str1);
            fprintf(fid,str3);
            fprintf(fid,str2);
        end
    end
    for i = 1:3;
        str = ['  ' num2str(76+i) '    1      1    0.00    0.00      0   0.30  0.30  0.00  0.50        75  75  75 x' '\r\n'];
        fprintf(fid,str);
    end
    
    
    for i = 1:numimages;
        if 79+i < 100
            str = ['  ' num2str(79+i) '    8           0.00    0.00      0                                 75  75  75 x   C:\\' set '\\' num2str(i) '.bmp' '\r\n'];
        else
            str = [' ' num2str(79+i) '    8           0.00    0.00      0                                 75  75  75 x   C:\\' set '\\' num2str(i) '.bmp' '\r\n'];
        end
        fprintf(fid,str);
    end
    fclose(fid);
    
    
    load('imagesequences1.mat'); %variable sequence matrix of image order
    fid = fopen([set '.CND'],'w+');
    
    line1 = ['COND#' ' ' 'BCKGND' ' ' 'TIMING' ' ' 'FIX_ID' ' ' '---COLOR-PALETTE---' ' ' 'TEST0' ' ' 'TEST1' '\r\n'];
    fprintf(fid,line1);
    
    clrchng = [];
    for i = 1:10;
        clrchng = [clrchng 1:ngrid*ngrid];
    end
    clrchng = clrchng(randperm(ngrid*ngrid*10));
    
    %--add 2 extra calibration trials at the beginging to get offset correct(to be ignored)---%
    background = num2str(-2);
    timing = num2str(2);
    fixationid = num2str(38);
    test0 = num2str(39);
    test1 = num2str(40);
    str = ['    ' num2str(1) '  ' background '      ' timing '      ' fixationid '                         ' test0 '     ' test1 '\r\n'];
    fprintf(fid,str);
    str = ['    ' num2str(2) '  ' background '      ' timing '      ' fixationid '                         ' test0 '     ' test1 '\r\n'];
    fprintf(fid,str);
    
    for iter = 1:50
        %calibratin trails%
        background = num2str(-2);
        timing = num2str(2);
        for i = 1:5
            clrchngtrial = clrchng((iter-1)*5+i);
            fixationid = num2str(1+(clrchngtrial-1)*3+1);
            test0 = num2str(1+(clrchngtrial-1)*3+2);
            test1 = num2str(1+(clrchngtrial-1)*3+3);
            
            if (iter-1)*6+i+2 < 10
                space1 = '    ';
            elseif (iter-1)*6+i+2 < 100
                space1 = '   ';
            else
                space1 = '  ';
            end
            if 1+(clrchngtrial-1)*3+1< 10
                space2 = '       ';
            else
                space2 = '      ';
            end
            if 1+(clrchngtrial-1)*3+2 < 10
                space3 = '                          ';
            else
                space3 = '                         ';
            end
            if 1+(clrchngtrial-1)*3+3 < 10
                space4 = '      ';
            else
                space4 = '     ';
            end
            
            str = [space1 num2str((iter-1)*6+i+2) '  ' background '      ' timing space2 fixationid space3 test0 space4 test1 '\r\n'];
            fprintf(fid,str);
        end
        
        %image trial%
        background = num2str(-4);
        timing = num2str(1);
        fixationid = num2str(-3);
        test0 = num2str(79+sequences(setnum,iter));
        if iter*6+2 < 10
            space1 = '    ';
        elseif iter*6+2 < 100
            space1 = '   ';
        else
            space1 = '  ';
        end
        str = [space1 num2str(iter*6+2) '  ' background '      ' timing '      ' fixationid '                         ' test0 '\r\n'];
        fprintf(fid,str);
    end
    fclose(fid);
end
%% [3] Generate Item files for wide screen setup(s)
% currently for  MP's setup. 
% CND files are the same. The only thing that changes is the background
% values in the first 7 lines.
numimages = 25;
for ssss = 2:5
    setnum = ssss;
    if setnum < 10
        set = ['WSMT00' num2str(setnum)];
    else
        set = ['WSMT0' num2str(setnum)];
    end
    
    fid = fopen([set '.itm'],'w+');
    
    spacingx = {'-12.0','-6.00',' 0.00',' 6.00',' 12.0'};
    spacingy = {'-8.00','-4.00',' 0.00',' 4.00',' 8.00'};
    ngrid = length(spacingx);
    
    line1 ='ITEM TYPE FILLED CENTERX CENTERY BITPAN HEIGHT WIDTH ANGLE OUTER INNER -R- -G- -B- C A ------FILENAME------';
    line2 ='  -4    1      1    0.00    0.00      0   1.00  1.00  0.00               0   0   0 x';
    line3 ='  -3   14      1    0.00    0.00      0   0.50  0.50  0.00             255 255 255 x';
    line4 ='  -2    1      1    0.00    0.00      0   0.20  0.20  0.00  0.50         0   0   0 x';
    line5 ='  -1    1      0    0.00    0.00      0   1.00  1.00  0.00             255 255 255 x';
    line6 ='   0    1      1    0.00    0.00      0   0.15  0.15  0.00             150 150 150 x';
    line7 ='   1    1      1    0.00    0.00      0   0.30  0.30  0.00  0.50       150 150 150 x';
    

    for line = 1:7
        fprintf(fid,[eval(['line' num2str(line)]) '\r\n']);
    end
    for i = 1:length(spacingx)
        for ii = 1:length(spacingy)
            num1 = num2str(1+3*ngrid*(i-1)+2*(ii-1)+ii);
            num3 = num2str(1+3*ngrid*(i-1)+2*(ii-1)+ii+1);
            num2 = num2str(1+3*ngrid*(i-1)+2*(ii-1)+ii+2);
            if str2double(num1) < 10
                header1 = ['   ' num1];
            else
                header1 = ['  ' num1];
            end
            if str2double(num2) < 10
                header2 = ['   ' num2];
            else
                header2 = ['  ' num2];
            end
            if str2double(num3) < 10
                header3 = ['   ' num3];
            else
                header3 = ['  ' num3];
            end
            str1 = [header1 '    1      1   ' spacingx{i} '   ' spacingy{ii} '      0   0.30  0.30  0.00  0.50       150 150 150 x' '\r\n'];
            str3 = [header3 '    1      1   ' spacingx{i} '   ' spacingy{ii} '      0   0.30  0.30  0.00  0.50       150 150 150 x' '\r\n'];
            str2 = [header2 '    1      1   ' spacingx{i} '   ' spacingy{ii} '      0   0.30  0.30  0.00  0.50       175 175 130 x' '\r\n'];
            fprintf(fid,str1);
            fprintf(fid,str3);
            fprintf(fid,str2);
        end
    end
    for i = 1:3;
        str = ['  ' num2str(76+i) '    1      1    0.00    0.00      0   0.30  0.30  0.00  0.50        75  75  75 x' '\r\n'];
        fprintf(fid,str);
    end
    
    
    for i = 1:numimages;
        if 79+i < 100
            str = ['  ' num2str(79+i) '    8           0.00    0.00      0                                 75  75  75 x   C:\\' set(2:end) '\\' num2str(i) '.bmp' '\r\n'];
        else
            str = [' ' num2str(79+i) '    8           0.00    0.00      0                                 75  75  75 x   C:\\' set(2:end) '\\' num2str(i) '.bmp' '\r\n'];
        end
        fprintf(fid,str);
    end
    fclose(fid);
end