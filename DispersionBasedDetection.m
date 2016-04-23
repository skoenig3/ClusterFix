%% I-MST: mininimum spanning tree method
tic

window_length =                   100; %50, window size

scm_image_dir = 'C:\Users\seth.koenig\Documents\MATLAB\BCRW Salience Model\SCM Image Sets\';
image_sets = {'Set006','Set007','Set008','Set009',...
    'SetE001','SetE002','SetE003','SetE004'};

samprate = 5/1000;
fltord = 60;
lowpasfrq = 30;
nyqfrq = 1000 ./ 2;
flt = fir2(fltord,[0,lowpasfrq./nyqfrq,lowpasfrq./nyqfrq,1],[1,1,0,0]);
buffer = 100/samprate/1000; %100 ms buffer for filtering 

for idt = [NaN 0.5 0.75 1 1.25];
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
                distances = [];
                x = eyedat{cndlop}(1,:)*24+400; %converts data from [-400,400] to [0,800]
                y = eyedat{cndlop}(2,:)*24+300; %converts data from [-300,300] to [0,600]
                x = [x(buffer:-1:1) x x(end:-1:end-buffer)]; %add buffer for filtering
                y = [y(buffer:-1:1) y y(end:-1:end-buffer)];   %add buffer for filtering
                x = resample(x,samprate*1000,1);%up sample to 1000 Hz
                y = resample(y,samprate*1000,1);%up sample to 1000 Hz
                x = filtfilt(flt,1,x);
                y = filtfilt(flt,1,y);
                x = x(101:end); %remove buffer for filtering but leave 100 ms for window buffer
                y = y(101:end); %remove buffer for filtering but leave 100 ms for window buffer
                
                fixationindexes = [];
                saccadeindexes = [];
                
                windows = round(length(x)/window_length); %insures at least a window of 100 ms
                for  ws = 1:windows;
                    % Computing Euclidean distance between each eye position inside the window
                    window_ind = window_length*(ws-1)+1:ws*window_length+1;
                    window_ind(window_ind > length(x)) = [];
                    window_len = length(window_ind);
                    X = zeros(window_len,2);
                    X(:,1) = x(window_ind);
                    X(:,2) = y(window_ind);
                    distance_matrix = squareform(pdist(X,'euclidean'));
                    distance_matrix(logical(eye(size(distance_matrix)))) = NaN;
                    vertices = zeros(window_len,1);    % List of visited vertices
                    distances = zeros(1,window_len-1);
                    vertices(1) = 1;                        % Mark first vertex as visited
                    count_v = 1;                            % Total amount of visited vertices
                    while (count_v<window_len)           % While some of the vertices is unvisited
                        min_distance = inf;                 % Searching for minimal distance betwenn visited and unvisited vertices
                        for i=1:window_len
                            if (vertices(i) == 1)
                                for j=1:window_len
                                    if (vertices(j) == 0)
                                        if( min_distance > distance_matrix( i,j ) )
                                            min_distance = distance_matrix( i,j );
                                            v1 = i;
                                            v2 = j;
                                        end
                                    end
                                end
                            end
                        end
                        vertices(v2) = 1;                   % Mark selected unvisited vertex as visited
                        distances(count_v) = min_distance;
                        count_v = count_v + 1;
                    end
                    if isnan(idt)
                        distthresh = mean(distances)+std(distances);
                    else
                        distthresh = idt;
                    end
                    fixationindexes = [fixationindexes window_ind(find(distances <= distthresh))];
                    saccadeindexes = [saccadeindexes window_ind(find(distances > distthresh))];
                end
                fixationindexes(fixationindexes > length(x)-100) = [];
                saccadeindexes(saccadeindexes > length(x)-100) = [];
               
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
                if ~isempty((sacind))
                    sacations = zeros(2,size(sacind,1));
                    for index=1:size(sacind,1)
                        rowsacind = sacind(index,:);
                        rowsacind(rowsacind == 0) = [];
                        saccadetimes(:,index) = [rowsacind(1);rowsacind(end)];
                    end
                else
                    saccadetimes = [];
                end
                
                %---Consolidate using duration thresholds---%
                tooshort = find(diff(saccadetimes,1) < 10); %10 ms duration threshold for saccades
                notbehav = [];
                for ii = 1:length(tooshort);
                    notbehav = [notbehav saccadetimes(1,tooshort(ii)):saccadetimes(2,tooshort(ii))];
                end
                fixationindexes = sort([fixationindexes notbehav]);
                
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
                for index=1:size(fixind,1)
                    rowfixind = fixind(index,:);
                    rowfixind(rowfixind == 0) = [];
                    fixationtimes(:,index) = [rowfixind(1);rowfixind(end)];
                end
                
                fixationtimes(:,(diff(fixationtimes,1) < 25))= []; %25 ms duration threshold
                fixationindexes = [];
                for ii = 1:size(fixationtimes,2);
                    fixationindexes = [fixationindexes fixationtimes(1,ii):fixationtimes(2,ii)];
                end
                
                saccadeindexes = 1:length(x)-100;
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
                if ~isempty((sacind))
                    sacations = zeros(2,size(sacind,1));
                    for index=1:size(sacind,1)
                        rowsacind = sacind(index,:);
                        rowsacind(rowsacind == 0) = [];
                        saccadetimes(:,index) = [rowsacind(1);rowsacind(end)];
                    end
                else
                    saccadetimes = [];
                end
                
                %---Return indexes to previous sampling rate & Calculate mean fixation position---%
                fixationtimes = fixationtimes;
                round5 = rem(fixationtimes,samprate*1000);
                round5(1,round5(1,:) > 0) = samprate*1000-round5(1,round5(1,:) > 0);
                round5(2,:) = - round5(2,:);
                fixationtimes = round((fixationtimes+round5)/5);
                fixationtimes(fixationtimes < 1) = 1;
                
                if ~isempty(saccadetimes)
                    saccadetimes = saccadetimes;
                    round5 = rem(saccadetimes,samprate*1000);
                    round5(1,:) = - round5(1,:);
                    round5(2,round5(2,:) > 0) = samprate*1000-round5(2,round5(2,:) > 0);
                    saccadetimes = round((saccadetimes+round5)/5);
                    saccadetimes(saccadetimes < 1) = 1;
                end
                
                x = eyedat{cndlop}(1,:)*24+400;
                y = eyedat{cndlop}(2,:)*24+300;
                saccadetimes(saccadetimes > length(x)) = length(x);
                fixationtimes(fixationtimes > length(x)) = length(x);
                
                fixations = zeros(size(fixationtimes));
                for ii = 1:size(fixationtimes,2)
                    fixations(1,ii) = mean(x(fixationtimes(1,ii):fixationtimes(2,ii)));
                    fixations(2,ii) = mean(y(fixationtimes(1,ii):fixationtimes(2,ii)));
                end
                
                fixationstats{cndlop}.fixationtimes = fixationtimes;
                fixationstats{cndlop}.fixations = fixations;
                fixationstats{cndlop}.saccadetimes = saccadetimes;
                fixationstats{cndlop}.XY = [x;y];
                
            end
            save([cortexfile(1:end-2) '_' cortexfile(end) '-MST-' num2str(idt) '.mat'],'fixationstats')
        end
    end
end
toc
%% I - DT: dispersion threshold method
tic
scm_image_dir = 'C:\Users\seth.koenig\Documents\MATLAB\BCRW Salience Model\SCM Image Sets\';
image_sets = {'Set006','Set007','Set008','Set009',...
    'SetE001','SetE002','SetE003','SetE004'};

IDTminduration = 25;
samprate = 5/1000;
fltord = 60;
lowpasfrq = 30;
nyqfrq = 1000 ./ 2;
flt = fir2(fltord,[0,lowpasfrq./nyqfrq,lowpasfrq./nyqfrq,1],[1,1,0,0]);
buffer = 100/samprate/1000; %100 ms buffer for filtering 

for idt = [0.5 0.75 1 1.25];
    
    IDTthresh = idt*24;
    
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
                x = eyedat{cndlop}(1,:)*24+400; %converts data from [-400,400] to [0,800]
                y = eyedat{cndlop}(2,:)*24+300; %converts data from [-300,300] to [0,600]
                x = [x(buffer:-1:1) x x(end:-1:end-buffer)]; %add buffer for filtering
                y = [y(buffer:-1:1) y y(end:-1:end-buffer)];   %add buffer for filtering
                x = resample(x,samprate*1000,1);%up sample to 1000 Hz
                y = resample(y,samprate*1000,1);%up sample to 1000 Hz
                x = filtfilt(flt,1,x);
                y = filtfilt(flt,1,y);
                x= x(101:end-75); %remove buffer for filtering but leave 25 ms for window buffer
                y= y(101:end-75); %remove buffer for filtering but leave 25 ms for window buffer

                fixationindexes = zeros(1,length(x));
                window_start = 1;
                window_end = IDTminduration;
                while window_end <= length(x);
                    IDT_max_x = max(x(window_start:window_end));
                    IDT_max_y = max(y(window_start:window_end));
                    IDT_min_x = min(x(window_start:window_end));
                    IDT_min_y = min(y(window_start:window_end));
                    result = abs( IDT_max_x - IDT_min_x ) + abs( IDT_max_y - IDT_min_y );
                    if result <= IDTthresh
                        if window_end-window_start + 1 == 25;
                            fixationindexes(window_start:window_end) = 1;
                        end
                        fixationindexes(window_end) = 1;
                        window_end = window_end +1;
                    else
                        window_start = window_start + 1;
                        window_end = window_start + IDTminduration - 1;
                    end
                end
                fixationindexes = find(fixationindexes);
                saccadeindexes = 1:length(x)-25;
                [~ , ~, ib] = intersect(fixationindexes,saccadeindexes);
                saccadeindexes(ib) = [];
                
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
                for index=1:size(fixind,1)
                    rowfixind = fixind(index,:);
                    rowfixind(rowfixind == 0) = [];
                    fixationtimes(:,index) = [rowfixind(1);rowfixind(end)];
                end
                
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
                if ~isempty((sacind))
                    sacations = zeros(2,size(sacind,1));
                    for index=1:size(sacind,1)
                        rowsacind = sacind(index,:);
                        rowsacind(rowsacind == 0) = [];
                        saccadetimes(:,index) = [rowsacind(1);rowsacind(end)];
                    end
                else
                    saccadetimes = [];
                end
                
                fixationtimes = fixationtimes;
                round5 = rem(fixationtimes,samprate*1000);
                round5(1,round5(1,:) > 0) = samprate*1000-round5(1,round5(1,:) > 0);
                round5(2,:) = - round5(2,:);
                fixationtimes = round((fixationtimes+round5)/5);
                fixationtimes(fixationtimes < 1) = 1;
                
                if ~isempty(saccadetimes)
                    saccadetimes = saccadetimes;
                    round5 = rem(saccadetimes,samprate*1000);
                    round5(1,:) = - round5(1,:);
                    round5(2,round5(2,:) > 0) = samprate*1000-round5(2,round5(2,:) > 0);
                    saccadetimes = round((saccadetimes+round5)/5);
                    saccadetimes(saccadetimes < 1) = 1;
                end
                
                x = eyedat{cndlop}(1,:)*24+400;
                y = eyedat{cndlop}(2,:)*24+300;
                saccadetimes(saccadetimes > length(x)) = length(x);
                fixationtimes(fixationtimes > length(x)) = length(x);
                
                fixations = zeros(size(fixationtimes));
                for ii = 1:size(fixationtimes,2)
                    fixations(1,ii) = mean(x(fixationtimes(1,ii):fixationtimes(2,ii)));
                    fixations(2,ii) = mean(y(fixationtimes(1,ii):fixationtimes(2,ii)));
                end
                
                fixationstats{cndlop}.fixationtimes = fixationtimes;
                fixationstats{cndlop}.fixations = fixations;
                fixationstats{cndlop}.saccadetimes = saccadetimes;
                fixationstats{cndlop}.XY = [x;y];
            end
            save([cortexfile(1:end-2) '_' cortexfile(end) '-IDT-' num2str(IDTthresh) '.mat'],'fixationstats')
        end
    end
end
toc