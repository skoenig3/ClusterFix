% Alternative Frequency Cluster_Fix

load('MEM_Pt011_block1_eyedat.mat')
eyedat = MEM_Pt011_block1_eyedat(:,2:3);
eyedat = [eyedat(2000:3000,1)';eyedat(2000:3000,2)'];
eyedat = {eyedat};

[fixationstats] = ClusterFixation_MJ(eyedat,1/125); 

%%
figure
for fi = 1:length(fixationstats)
x = fixationstats{fi}.XY(1,:);
y = fixationstats{fi}.XY(2,:);


plot(x(25:end-25),y(25:end-25),'g') %exclude first and last 200 ms
hold on
for i = 1:size(fixationstats{fi}.fixationtimes,2);
    times = fixationstats{fi}.fixationtimes(1,i):fixationstats{fi}.fixationtimes(2,i);
    plot(x(times),y(times),'r')
    fix = fixationstats{fi}.fixations(:,i);
    plot(fix(1),fix(2),'k*')
end
legend('Saccades','Fixations')
hold off
pause(3)
end