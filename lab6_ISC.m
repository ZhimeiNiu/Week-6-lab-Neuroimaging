clear all

%lab 6 ISC
%load file (create toy dataset)
n = 20;
r = 5;
t = 50;
s = 2;
behavior = rand([n,n]);
roi_data = rand([n,r,t,s]);
%% Script settings:
temporal = 1;
dynamic = 1;
spatial = 1;
intra_sub = 1;
ISFC = 1;
behave_isc = 1;
%% 1.Temporal ISC

%% do pairwise temporal ISC
if temporal

pairwise_temporal_ISC = zeros([n,n,r]);
%average across timeseries
mean_tseries = mean(roi_data,4);

%do correlation
for i=1:n
    for j=1:n
        for x=1:r
            pairwise_temporal_ISC(i,j,x)= corr(squeeze(mean_tseries(i,x,:)),squeeze(mean_tseries(j,x,:)));
        end
    end
end

%% do leave-one-out ISC
loo_temporal_ISC =zeros([n,r]);

for i=1:n
    for x=1:r
        %leave one out (the sub itself)
        mean_tseries2 = mean_tseries
        mean_tseries2(i,:,:)=[]; 
        loo_temporal_ISC(i,x)= corr(squeeze(mean_tseries(i,x,:)), squeeze(mean(mean_tseries2(:,x,:))));
    end
end
end
%% Dynamic ISC
if dynamic
%window size and step size
wins = 10;
step = 10;
w = (t-wins)/step+1
%create dynamic_ISC
loo_dynamic_ISC= zeros([n,r,w]);

for i=1:n
    for x=1:r
        for y = 1:w
        %leave one out
        mean_tseries3 = mean_tseries
        mean_tseries3(i,:,:)=[];
        loo_dynamic_ISC(i,x,y)= corr(squeeze(mean_tseries(i,x,(y-1)*w+1:y*w)),squeeze(mean(mean_tseries3(:,x,(y-1)*w+1:y*w))));
        end
    end
end

figure
plot(10:10:50,squeeze(loo_dynamic_ISC(4,2,:)))
xlabel('Timecourse(s)')
ylabel('Dynamic ISC')

end
%% (3) Spatial ISC
if spatial

loo_spatial_ISC = zeros([n,1]);

for i=1:n
        mean_tseries4  = mean_tseries;
        mean_tseries4(i,:,:)=[]; 
        loo_spatial_ISC(i)= corr(mean(mean_tseries(i,:,:),3)',mean(mean(mean_tseries4(:,:,:),3))') ;
end
end
%% Intra-subject correlation
if intra_sub
intrasubject_temporal_ISC = zeros([n,r]);
for i=1:n
    for j=1:r
        intrasubject_temporal_ISC(i,j)= corr(squeeze(roi_data(i,j,:,1)),squeeze(roi_data(i,j,:,2))) ;
    end
end
end
%% Inter-subject functional connectivity
% leave one out ISC
if ISFC
loo_ISFC =zeros([n,r,r]);
for i=1:n
    for x=1:r
        for y=1:r
            mean_tseries5  = mean_tseries;
            mean_tseries5(i,:,:)=[]; 
            loo_ISFC(i,x,y)= corr(squeeze(mean_tseries(i,x,:)), squeeze(mean(mean_tseries5(:,y,:))));
        end
    end
end
end
%% Inter-subject similarity to behavior
if behave_isc
behavioral_similarity = zeros([n,n]);

for x=1:n
    for y=1:n
        behavioral_similarity(x,y) = sum(behavior(x,:)-behavior(y,:).^2);
    end
end

roi_sim = zeros(5,1);
pvalue = zeros(5,1);
%behavioral similarity convert matrix 
matrix_be = ones(size(behavioral_similarity,1));
matrix_half = tril(matrix_be,-1);
sum_vector = behavioral_similarity(:);
behavioral_similarity  = sum_vector(matrix_half(:)>0);

for i=1:r

    [roi_sim(i),pvalue(i)]=corr(behavioral_similarity,shape(squeeze(pairwise_temporal_ISC(:,:,i))));
end

end
%the hypothesis should be that the ISC and behavioral response are not
%correlated since the data are all random numbers

%I correlate the ISC and similarity between behavioral response for each
%roi

%result
% there is no correlation for the behavioral similarity to ISC
% p values are larger than 0.05

function [output] = shape(input)
matrix_be2 = ones(size(input,1));
matrix_half2 = tril(matrix_be2,-1);
sum_vector2 = input(:);
output  = sum_vector2(matrix_half2(:)>0);
end
