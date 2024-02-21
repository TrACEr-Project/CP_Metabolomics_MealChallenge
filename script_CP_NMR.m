% This script shows the workflow to fit a CP model to NMR measurements
% of plasma samples collected during a challenge test in the manuscript:
%
%   S. Yan, L. Li, D. Horner, P. Ebrahimi, B. Chawes, L. O. Dragsted, 
%   M. A. Rasmussen, A. K. Smilde,  E. Acar, Characterizing human postprandial 
%   metabolic response using multiway data analysis, 2023 
%   doi: https://doi.org/10.1101/2023.08.31.555521
%
% The script assumes that the data is stored as a dataset object
% (https://eigenvector.com/software/dataset-object/)
%
% Needed toolboxes:
% Tensor toolbox (https://www.tensortoolbox.org/)
% Poblano toolbox (https://github.com/sandialabs/poblano_toolbox)

%% load data and prepare the data for the analysis
load('data.mat')
% NMR is of size 299 subjects by 8 time points by 252 features (i.e., Nightingale
% measurements + insulin + cpeptide)
exc_id = [69 75]; % remove Phe and Glycerol due to missing data
X      = NMR(:,:,setdiff(1:size(NMR,3), exc_id));

% Select the 161 features (metabolites + insulin +c-peptide) included in the study
id = find(X.class{3,1}>0);
X  = X(:,:,id);

% select males (id_gender = 2) or females (id_gender = 1)
id_gender = find(X.class{1,1}==2);
YY        = X(id_gender,:,:); %1:females, 2:males
Meta      = Metainfo(id_gender,:); % additional information available about the subjects

% remove outliers
X = YY;
outlier  = [142 79 342 312]; %male
%outlier = [90 335 250]; %female
[~,id,~] = intersect(str2num(X.label{1,1}), outlier);
inc      = setdiff(1:size(X,1),id);
Xfinal   = X(inc,:,:);
Meta     = Meta(inc,:);

% merge BMI groups into two groups
cc = Xfinal.class{1,2};
cc(find(Xfinal.class{1,2}==2 | Xfinal.class{1,2}==3))=2; %obese and overweight
cc(find(Xfinal.class{1,2}==1 | Xfinal.class{1,2}==4))=1; %underweight and normal
Xfinal.class{1,11}=cc;
Xfinal.classname{1,11}='BMI 2-class';

% fasting and T0-corrected data 
Xfast  = squeeze(Xfinal(:,1,:));
XT0    = take_diff(Xfinal); 

%% fit CP model to T0-corrected using cp_wopt
Xpre = preprocess_centerscale(XT0,1,1);

W  = tensor(~isnan(Xpre.data));
XX = tensor(Xpre.data);
XX(find(W==0))=0;

% optimization parameters 
options = ncg('defaults');
options.MaxFuncEvals = 10000;
options.MaxIters     = 10000;
options.StopTol      = 1e-10;
options.RelFuncTol   = 1e-8;
options.DisplayIters = 100;
nb_start = 10;
R = 2;
for i=1:nb_start
    [Fac{i},FacInit{i},out{i}] = cp_wopt(XX,W, R, 'init','randn','opt','ncg','opt_options',options);
    f(i) = out{i}.f;
end
[ff, index] = sort(f,'ascend');
if out{index(1)}.ExitFlag==3 || out{index(1)}.ExitFlag==0 % algorithm stops due to relative change in function value or the gradient condition
    flag_stop = 1;
    Fac_cp = Fac{index(1)};
    out_cp = out{index(1)};
else 
    flag_stop = 0;
end
%% Replicability for different number of components
R = 4;
for r = 1:R
    allscores(:,r) = CP_replicability(XT0,r);
end
plot_allscores(allscores)

            
%% Example showing how to check the results, e.g., looking at subject scores, t-test, correation with variables of interest
S      = Fac_cp.U{1}; % Subject scores
[h, p] = ttest2(S(Xfinal.class{1,11}==1,:),S(Xfinal.class{1,11}==2,:),'alpha', 0.05, 'vartype','unequal');

ind = find(h==1);
for i = 1:length(ind)
    x{1}   = S(Xfinal.class{1,11}==1,ind(i)); x{2}=S(Xfinal.class{1,11}==2,ind(i));
    xx{i}  = [x{1};x{2}];
    g      = [ones(size(x{1},1),1);2*ones(size(x{2},1),1)];
    boxplot(xx{i}, g); title(strcat('Large Dataset -comp', num2str(ind(i))));
    set(gca, 'XTick',1:1:2, 'XTickLabel',{'Lower BMI', 'Higher BMI'})
    clear x
end

% How do the scores correlate with variables of interest?
comp_id = 2; %relevant component id
figure
for i=1:size(Meta,2)
    idnan = isnan(Meta.data(:,i));
    if sum(idnan)<50
        [C(i), p(i)]=corr(S(idnan==0,comp_id),Meta.data(find(idnan==0),i));
    else
        C(i)= NaN;
        p(i)= NaN;
     end
end
