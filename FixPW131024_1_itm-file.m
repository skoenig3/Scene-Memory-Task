%% SMT001 was accidentally run on an inccorectly compiled smt.sav file with wrong 
% CND file so have to "remake" file to match color change trials. Picture
% trials were ok

% First run getSMTeyedat unitl %Test for errors% in debug mode
% insert following lines after %---For Calibration with Eye tracking data with cp2tform---%
oldcnd = cnd-1000;
newcnd = zeros(length(oldcnd),2);

% insert following lines at the end of the first for loop after %---For Calibration with Eye tracking data with cp2tform---%

mx = mean(eog_arr(oddind,per(k).event));
my = mean(eog_arr(evenind,per(k).event));
if mx < -1000
    newcnd(k,1) = 1;
elseif mx < -500
    newcnd(k,1) = 2;
elseif mx < 500
    newcnd(k,1) = 3;
elseif mx < 1000
    newcnd(k,1) = 4;
else
    newcnd(k,1) = 5;
end
if my < -1000
    newcnd(k,2) = 1;
elseif my < -500
    newcnd(k,2) = 2;
elseif my < 400
    newcnd(k,2) = 3;
elseif my < 1000
    newcnd(k,2) = 4;
else
    newcnd(k,2) = 5;
end
%% replace item numbers in PWSMT001.CND

for k = 1:size(newcnd,1);
    for line = 2:size(itmfil,1);
        C = textscan(itmfil(line,:),'%d');
        if C{1}(4) == spacex(newcnd(k,1)) & C{1}(5) == spacey(newcnd(k,2))
            break
        end
    end
    background = num2str(-2);
    timing = num2str(2);
    fixationid = num2str(line-6);
    test0 = num2str(line-6+1);
    test1 = num2str(line-6+2);
    
    if  oldcnd(k) < 10
        space1 = '    ';
    elseif oldcnd(k) < 100
        space1 = '   ';
    else
        space1 = '  ';
    end
    if fixationid < 10
        space2 = '        ';
    else
        space2 = '    ';
    end
    if test0 < 10
        space3 = '                         ';
    else
        space3 = '                        ';
    end
    if test1 < 10
        space4 = '     ';
    else
        space4 = '    ';
    end
    
    str = [space1 num2str(oldcnd(k)) '  ' background '       ' timing space2 fixationid space3 test0 space4 test1];
    
    for i = 1:58-length(str)
        str = [str ' '];
    end
    
    cndfil(oldcnd(k)+1,:) = str;
end

