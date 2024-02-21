function Xpre = preprocess_centerscale(X, flag_scale, flag_center)
% Given the tensor X: subjects by time by metabolites 
%center across subject, scale within metabolites mode using rms

% preprocessing
Xpre = X;
if flag_center
    temp           = X.data(:,:);
    temp_centered  = temp - repmat(nanmean(temp),size(temp,1),1);
    Xpre.data      = reshape(temp_centered, size(X));
end
if flag_scale && length(size(Xpre))==3
    for j=1:size(Xpre,3)
        temp = squeeze(Xpre.data(:,:,j));
%        temp = temp/nanstd(temp(:)); %scaling with std
        rms  = sqrt(nanmean(temp(:).^2));
        temp = temp/rms;
        XX(:,:,j)= temp;
    end
    Xpre.data = XX;
elseif flag_scale && length(size(Xpre))==2
    for j=1:size(Xpre,2)
        temp = squeeze(Xpre.data(:,j));
%        temp = temp/nanstd(temp(:)); %scaling with std
        rms  = sqrt(nanmean(temp(:).^2));
        temp = temp/rms;
        XX(:,j)= temp;
    end
     Xpre.data = XX;
end


