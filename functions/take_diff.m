function Xdiff = take_diff(X)

XX = X;
for i=2:size(X,2)
    XX.data(:,i-1,:) = X.data(:,i,:)-X.data(:,1,:);
end
Xdiff = XX(:,1:end-1,:);
Xdiff.axisscale{2} = XX.axisscale{2}(2:end);
Xdiff.label{2}     = XX.label{2}(2:end,:);
