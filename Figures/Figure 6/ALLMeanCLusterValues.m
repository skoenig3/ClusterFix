% Mean Cluster Values
scm_image_dir = 'C:\Users\seth.koenig\Documents\MATLAB\BCRW Salience Model\SCM Image Sets\';
image_sets = {'Set006','Set007','Set008','Set009',...
    'SetE001','SetE002','SetE003','SetE004'};
tags = {'MP','TT','JN','IW'};
imageX = 800; imageY = 600;
fixclusters = cell(1,length(tags));
sacclusters = cell(1,length(tags));
eyedatafiles = NaN(1,length(tags));
for SET = 1:length(image_sets);
    SETNUM = image_sets{SET};
    cd([scm_image_dir SETNUM])
    matfiles = what;
    for i = 1:length(matfiles.mat);
        str = strfind(matfiles.mat{i},'fixation');
        if ~isempty(str)
            if isempty(strfind(matfiles.mat{i},'thresh'))
                for ii =1:length(tags);
                    if strfind(matfiles.mat{i},tags{ii})
                        eyedatafiles(ii) = i;
                    end
                end
            end
        end
    end
    for ii = 1:length(tags);
        if ~isnan(eyedatafiles(ii))
            load(matfiles.mat{eyedatafiles(ii)});
            for cndlop = 1:2:length(fixationstats)
                fixclusters{ii} = [fixclusters{ii} ...
                    fixationstats{cndlop}.MeanClusterValues(1,2:2:end)'];
                sacclusters{ii} = [sacclusters{ii} ...
                    fixationstats{cndlop}.MeanClusterValues(2,2:2:end)'];
            end
        end
    end
end
%% Plot angle, max velocity and distance-do not use without chaning lines 31 and 33 above
for i = 1:length(tags)
    figure
    hold on
    for ii = 1:size(fixclusters{i},2)
        plot3(fixclusters{i}(1,ii),fixclusters{i}(2,ii)/24,fixclusters{i}(3,ii),'r.')
    end
    for ii = 1:size(sacclusters{i},2)
        plot3(sacclusters{i}(1,ii),sacclusters{i}(2,ii)/24,sacclusters{i}(3,ii),'g.')
    end
    view(3)
    ylabel('Average Max Velocity (dva/sec)')
    xlabel('Average Mean Distance (dva)')
    zlabel('Average |Angle| (degrees)')
end
%% plot max acceleration, mean velcoity, and rotation
for i = 1:length(tags)
    figure
    hold on
    for ii = 1:size(fixclusters{i},2)
        plot3(fixclusters{i}(1,ii)/24,fixclusters{i}(2,ii)/24,200*fixclusters{i}(3,ii),'r.')
    end
    for ii = 1:size(sacclusters{i},2)
        plot3(sacclusters{i}(1,ii)/24,sacclusters{i}(2,ii)/24,200*sacclusters{i}(3,ii),'g.')
    end
    view(3)
    xlabel('Average Max Accleration (dva/sec^{2})')
    ylabel('Average Vel (dva/sec)')
    zlabel('Average Rotation (degrees/sample)')
end
%% Cluster Clustered Data
bothclusters = [fixclusters{4} sacclusters{4}];
T = kmeans(bothclusters',2,'replicate',5);

clr = ['rg'];
figure
hold on
for i = 1:length(T);
    plot3(bothclusters(1,i),bothclusters(2,i),bothclusters(3,i),[clr(T(i)) '.'])
end
hold off

% Support Vector Machine (SVM) classifier
correct = NaN(length(tags),100);
totalpoints = NaN(length(tags),100);
for t = 1:length(tags)
    for iter = 1:100;
        data = [fixclusters{t}'; sacclusters{t}'];
        groups = [true(1,size(fixclusters{t},2)) false(1,size(sacclusters{t},2))]';
        P = cvpartition(groups,'Holdout',0.75); % 0.75 uses 25% of points to train and 75% to test
        SVMstruct = svmtrain(data(P.training,:),groups(P.training));%,'showplot','true');
        testresults = svmclassify(SVMstruct,data(P.test,:));%,'showplot','true');
        correct(t,iter) = sum(groups(P.test) == testresults);
        totalpoints(t,iter) = sum(P.test);
    end
end
lowestaccuracy = 100*min(correct./totalpoints)
averageaccuracy = 100*mean(correct./totalpoints)
