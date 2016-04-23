% Boot strapping method for analyzing consistency of Cluster Fix across
% multiple applications of the algoirth to the same scan paths.

% Only computers Raw data. CF_Reliability_Analysis analyzes raw data.

% Can run for 100 iterations (numiters) but causes lots of crashing and can
% take 48+ hours to run. 100 iterations takes like 6 hours and doesn't
% experience as many issues. 

addpath('C:\Users\seth.koenig\Documents\MATLAB\ClusterFix\');

eyedat  = getEyeData('IW100709.1',.005);
numiters = 100;
allclassifications = cell(1,72);
for trial = 1:2:72;
    eyetrace = eyedat(trial);
    
    classification = zeros(numiters,size(eyetrace{1},2));
    for i = 1:numiters
        fixstats = ClusterFixation_Final(eyetrace);
        fixationtimes = fixstats{1}.fixationtimes;
        for f = 1:size(fixationtimes,2);
            classification(i,fixationtimes(1,f):fixationtimes(2,f)) = 1;
        end
    end 
    
    allclassifications{trial} = classification;
    
end
save('C:\Users\seth.koenig\Documents\MATLAB\ClusterFix\Figures\Supplemental Figure 4-Reliability of CF\ReliabilityData_IW','allclassifications')
%% Same Boot strapping method as above for measuring Consistency of Cluster Fix
% but now using k-means ++ algoirthm instead of standard k-means algorithm
% that comes with matlab.

addpath('C:\Users\seth.koenig\Documents\MATLAB\ClusterFix\');

eyedat  = getEyeData('MP100712.1',.005);
numiters = 100;
allclassifications = cell(1,72);
for trial = 1:2:72;
    eyetrace = eyedat(trial);
    
    classification = zeros(numiters,size(eyetrace{1},2));
    for i = 1:numiters
        fixstats = ClusterFixation_PP(eyetrace);
        fixationtimes = fixstats{1}.fixationtimes;
        for f = 1:size(fixationtimes,2);
            classification(i,fixationtimes(1,f):fixationtimes(2,f)) = 1;
        end
    end 
    
    allclassifications{trial} = classification;
    
end
save('C:\Users\seth.koenig\Documents\MATLAB\ClusterFix\Figures\Supplemental Figure 4-Reliability of CF\ReliabilityData_MP_PP','allclassifications')
%% Same Boot strapping method as above for measuring Consistency of Cluster Fix
% but now using changing the number of replications performed by k-means

addpath('C:\Users\seth.koenig\Documents\MATLAB\ClusterFix\');

eyedat  = getEyeData('MP100712.1',.005);
numiters = 100;
allclassifications = cell(1,72);
for trial = 1:2:72;
    eyetrace = eyedat(trial);
    
    classification = zeros(numiters,size(eyetrace{1},2));
    for i = 1:numiters
        fixstats = ClusterFixation_Rep(eyetrace,5/1000,25);
        fixationtimes = fixstats{1}.fixationtimes;
        for f = 1:size(fixationtimes,2);
            classification(i,fixationtimes(1,f):fixationtimes(2,f)) = 1;
        end
    end 
    
    allclassifications{trial} = classification;
    
end
save('C:\Users\seth.koenig\Documents\MATLAB\ClusterFix\Figures\Supplemental Figure 4-Reliability of CF\ReliabilityData_MP_25','allclassifications')