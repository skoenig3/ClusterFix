%% Code runs a velocity and acceleration based detection algorithm for saccades and fixations.
% Velocity threshold and acceleration threhsold default at mean + 1*std.
% Velocity threshold finds consecutive points above threshold. Each set of
% consecutive points must contain at least one time point about the
% acceleration threhsold. Local re-valuation finds start and end of
% saccades based on acceleration threshold by finding time points 50 ms before
%and after the saccade above the acceleration threshold. Following re-evaluation,
% a 10 ms saccade duration threshold then a 25 ms fixation duration threshold are applied.
% This method follows Cluster Fix[ation] as closely as possible using
% velocity and acceleration thresholds. Written January 2013 by Seth Koenig.
scm_image_dir ='C:\Users\seth.koenig\Documents\MATLAB\BCRW Salience Model\SCM Image Sets\';
image_sets = {'Set006','Set007','Set008','Set009',...
    'SetE001','SetE002','SetE003','SetE004'};

samprate = 5/1000;
fltord = 60;
lowpasfrq = 30;
nyqfrq = 1000 ./ 2;
flt = fir2(fltord,[0,lowpasfrq./nyqfrq,lowpasfrq./nyqfrq,1],[1,1,0,0]);
buffer = 100/samprate/1000; %100 ms buffer for filtering

velocity = [NaN 30 75 100]; %NaN is mean + std of an individual scan path
acceleration = [NaN 8000 5000 8500]; %NaN is mean + std

for values = 1:length(velocity)
    for imset = 1:length(image_sets);
        dirName = [scm_image_dir image_sets{imset}];
        cd(dirName)
        dirData = dir(dirName);
        dirIndex = [dirData.isdir];
        fileList = {dirData(~dirIndex).name}';
        eyeindexes = [];
        for i = 1:length(fileList)
            period = strfind(fileList{i},'.');
            if (length(fileList{i}) == period(end)+1) && ...
                    (double(fileList{i}(period+1)) <=57)
                eyeindexes = [eyeindexes i];
            end
        end
        for i = 1:length(eyeindexes)
            cortexfile = fileList{eyeindexes(i)};
            eyedat  = getEyeData(cortexfile,samprate);
            fixationstats = cell(1,length(eyedat));
            for cndlop = 1:length(eyedat)
                x = eyedat{cndlop}(1,:)*24+400; %converts dva to pixel and data from [-400,400] to [0,800]
                y = eyedat{cndlop}(2,:)*24+300; %converts dva to pixel and from [-300,300] to [0,600]
                x = [x(buffer:-1:1) x x(end:-1:end-buffer)]; %add buffer for filtering
                y = [y(buffer:-1:1) y y(end:-1:end-buffer)];   %add buffer for filtering
                x = resample(x,samprate*1000,1);%up sample to 1000 Hz
                y = resample(y,samprate*1000,1);%up sample to 1000 Hz
                xss = filtfilt(flt,1,x);
                yss = filtfilt(flt,1,y);
                x = x(101:end-100);
                y = y(101:end-100);
                xss = xss(101:end-100); %remove buffer after filtering
                yss = yss(101:end-100); %remove buffer after filtering
                
                velx = diff(xss);
                vely = diff(yss);
                vel = sqrt(velx.^2+vely.^2);
                accel = abs(diff(vel));
                vel = vel(1:end-1);
                
                if isnan(velocity(values))
                    velthresh = mean(vel)+1*std(vel);
                    accelthresh =mean(accel)+1*std(accel);
                else
                    velthresh = velocity(values)/1000*24;% first # in degrees per sec
                    accelthresh = acceleration(values)/1000/1000*24;%first # in degrees per sec per sec
                end
                
                sacbeg = find(diff(vel > velthresh) > 0);
                sacend = find(diff(vel > velthresh) < 0);
                
                if vel(end)>=velthresh
                    if vel(1)> velthresh
                        tempbeg =[1 sacbeg];
                        tempend =[sacend length(accel)];
                    else
                        tempbeg =sacbeg;
                        tempend =[sacend length(accel)];
                    end
                else
                    if vel(1)> velthresh
                        tempbeg =[1 sacbeg];
                        tempend =sacend;
                    else
                        tempbeg = sacbeg;
                        tempend = sacend;
                    end
                end
                notsacs = [];
                for ii = 1:length(tempbeg);
                    if ~any(accel(tempbeg(ii):tempend(ii)) > accelthresh)
                        notsacs = [notsacs ii];
                    end
                end
                tempend(notsacs) = []; tempbeg(notsacs) = [];
                for ii = 1:length(tempbeg);
                    ind = tempbeg(ii)-50:tempbeg(ii)-1;
                    ind(ind < 1) = [];
                    ind(ind > length(accel)) = [];
                    prioraccel = accel(ind);
                    extrapiorsac = prioraccel > accelthresh;
                    nolongersac = find(extrapiorsac == 0);
                    if ~isempty(nolongersac)
                        ind(1:nolongersac(end)) = [];
                        if length(ind) >= 1
                            tempbeg(ii) = ind(1);
                        end
                    elseif ~isempty(ind)
                        tempbeg(ii) = ind(1);
                    end
                    
                    ind = tempend(ii)+1:tempend(ii)+50;
                    ind(ind < 1) = [];
                    ind(ind > length(accel)) = [];
                    postaccel = accel(ind);
                    extrapostsac = postaccel > accelthresh;
                    nolongersac = find(extrapostsac == 0);
                    if ~isempty(nolongersac)
                        ind(nolongersac(1):end) = [];
                        if length(ind) >= 1
                            tempend(ii) = ind(end);
                        end
                    elseif ~isempty(ind)
                        tempend(ii) = ind(end);
                    end
                end
                notsacs = find(tempend-tempbeg < 10);
                tempend(notsacs) = []; tempbeg(notsacs) = [];
                saccadetimes = [tempbeg;tempend];
                fixationtimes = [];
                for ii = 1:size(saccadetimes,2);
                    if ii == 1;
                        if saccadetimes(1,ii) > 1;
                            fixationtimes(:,1) = [1;saccadetimes(1,ii)-1];
                        else
                            fixationtimes(:,1) = [saccadetimes(2,1)+1;saccadetimes(1,2)-1];
                        end
                    elseif ii == size(saccadetimes,2);
                        fixationtimes(:,ii) = [saccadetimes(2,ii-1)+1,saccadetimes(1,ii)-1];
                        if saccadetimes(2,ii) < length(accel);
                            fixationtimes(:,ii+1) = [saccadetimes(2,ii)+1,length(accel)];
                        end
                    else
                        fixationtimes(:,ii) = [saccadetimes(2,ii-1)+1,saccadetimes(1,ii)-1];
                    end
                end
                notfixes = find(diff(fixationtimes) < 25);
                fixationtimes(:,notfixes) = [];
                
                fixationindexes = [];
                for ii = 1:size(fixationtimes,2);
                    fixationindexes = [fixationindexes fixationtimes(1,ii):fixationtimes(2,ii)];
                end
                
                dind = diff(fixationindexes);
                gaps =find(dind > 1);
                fixind = zeros(length(gaps),50);
                if ~isempty(gaps)
                    for gapind = 1:length(gaps)+1;
                        if gapind == 1;
                            temp = fixationindexes(1:gaps(gapind));
                        elseif gapind == length(gaps)+1
                            temp = fixationindexes(gaps(gapind-1)+1:end);
                        else
                            temp = fixationindexes(gaps(gapind-1)+1:gaps(gapind));
                        end
                        fixind(gapind,1:length(temp)) = temp;
                    end
                else
                    fixind =  fixationindexes;
                end
                fixationtimes = zeros(2,size(fixind,1));
                fixations = zeros(2,size(fixind,1));
                for index=1:size(fixind,1)
                    rowfixind = fixind(index,:);
                    rowfixind(rowfixind == 0) = [];
                    fixationtimes(:,index) = [rowfixind(1);rowfixind(end)];
                    fixations(:,index) = [mean(x(rowfixind(1):rowfixind(end)));...
                        mean(y(rowfixind(1):rowfixind(end)))];
                end
                
                saccadeindexes = 1:length(accel);
                [~ , ~, ib] = intersect(fixationindexes,saccadeindexes);
                saccadeindexes(ib) = [];
                
                dind = diff(saccadeindexes);
                gaps =find(dind > 1);
                sacind = zeros(length(gaps),50);
                if ~isempty(gaps)
                    for gapind = 1:length(gaps)+1;
                        if gapind == 1;
                            temp = saccadeindexes(1:gaps(gapind));
                        elseif gapind == length(gaps)+1
                            temp = saccadeindexes(gaps(gapind-1)+1:end);
                        else
                            temp = saccadeindexes(gaps(gapind-1)+1:gaps(gapind));
                        end
                        sacind(gapind,1:length(temp)) = temp;
                    end
                else
                    sacind =  saccadeindexes;
                end
                saccadetimes = zeros(2,size(sacind,1));
                sacations = zeros(2,size(sacind,1));
                for index=1:size(sacind,1)
                    rowsacind = sacind(index,:);
                    rowsacind(rowsacind == 0) = [];
                    saccadetimes(:,index) = [rowsacind(1);rowsacind(end)];
                end
                
                round5 = rem(fixationtimes,samprate*1000);
                round5(1,round5(1,:) > 0) = samprate*1000-round5(1,round5(1,:) > 0);
                round5(2,:) = - round5(2,:);
                fixationtimes = round((fixationtimes+round5)/5);
                fixationtimes(fixationtimes < 1) = 1;
                
                round5 = rem(saccadetimes,samprate*1000);
                round5(1,:) = - round5(1,:);
                round5(2,round5(2,:) > 0) = samprate*1000-round5(2,round5(2,:) > 0);
                saccadetimes = round((saccadetimes+round5)/5);
                saccadetimes(saccadetimes < 1) = 1;
                
                x = eyedat{cndlop}(1,:)*24+400;
                y = eyedat{cndlop}(2,:)*24+300;
                saccadetimes(saccadetimes > length(x)) = length(x);
                fixationtimes(fixationtimes > length(x)) = length(x);
                
                fixationstats{cndlop}.fixationtimes = fixationtimes;
                fixationstats{cndlop}.fixations = fixations;
                fixationstats{cndlop}.saccadetimes = saccadetimes;
                fixationstats{cndlop}.XY = [x;y];
            end
            if isnan(velocity(values))
                save([cortexfile(1:end-2) '_' cortexfile(end) '-threshmn.mat'],'fixationstats')
            else
                save([cortexfile(1:end-2) '_' cortexfile(end) '-threshv' num2str(velocity(values)) ...
                    'a' num2str(acceleration(values)) '.mat'],'fixationstats')
            end
        end
    end
end