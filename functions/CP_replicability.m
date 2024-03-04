function allscores = CP_replicability(X,R)
% This function checks the replicability of the CP model by 
%   Step1: forming a random ten-fold of the samples in the first mode of data tensor X 
%   Step2: leave one fold of the data out 
%   Step3: Fit an R-component CP model to the remaining 90% - repeat this
%   by leaving out each fold once
%   Step4: Compute the factor match scores between the best runs of the CP models (in the second and third modes). 
%   Step5: Repeat Step1-4 ten times
% Return in total 450 FMS scores.

%% Reorder subjects
sub_low    = find(X.class{1,2}==1 | X.class{1,2}==4); % Lower BMI subjects
sub_high   = find(X.class{1,2}==2 | X.class{1,2}==3); % Higher BMI subjects
index_perm = [sub_low,sub_high];
X          = X(index_perm,:,:);

%% optimization parameters
options = ncg('defaults');
options.MaxFuncEvals = 10000;
options.MaxIters     = 10000;
options.StopTol      = 1e-10;
options.RelFuncTol   = 1e-8;
options.DisplayIters = 100;

%% form splits and fit CP
nb_splits  = 10;
repeats    = 10;
allscores  = [];
for r = 1:repeats
    index_perm = [randperm(length(sub_low)),randperm(length(sub_high))+length(sub_low)];
    for s = 1:nb_splits
        index_rem = index_perm(s:10:end); %randomly remove 1/10 subjects (equally many in each group)
         % data with remaining subjects
        XX{s} = X(setdiff(index_perm, index_rem),:,:);     
        % preprocess the data: centering and scaling and imputing missing values
        XXprep = preprocess_centerscale(XX{s},1,1);
        [Fac{s}, flags(r,s)] = fitmodel(XXprep, R, options);
    end
    fms{r}    = compute_pairwise(Fac);
    allscores = [allscores;fms{r}];
    clear index Fac XX
end
eval(strcat('save allscores', num2str(R),'.mat'))



function [Fac_cp, flag_stop] = fitmodel(Y, R, options)

if R>2
    nb_starts = 50;
else
    nb_starts = 20;
end

W  = tensor(~isnan(Y.data));
XX = tensor(Y.data);
XX(find(W==0))=0;


for i=1:nb_starts
    [Fac{i},FacInit{i},out{i}] = cp_wopt(XX,W, R, 'init','randn','opt','ncg','opt_options',options);
    f(i) = out{i}.f;
end
[ff, index] = sort(f,'ascend');
Fac_cp = Fac{index(1)};
out_cp = out{index(1)};
if out{index(1)}.ExitFlag==3 || out{index(1)}.ExitFlag==0 % algorithm stops due to relative change in function value or the gradient condition
   flag_stop = 1;
else
   flag_stop = 0;
end



% compute_pairwise
function fms = compute_pairwise(Fac)

nb_runs = length(Fac);
for i=1:nb_runs
    for j=i+1:nb_runs
        Sim(i,j)= score(ktensor(Fac{i}.U([2:3])), ktensor(Fac{j}.U([2:3])));
    end
end
fms = Sim(find(Sim>0));
