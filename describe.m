function summary = describe(X)
[rx,cx] = size(X);
summary = struct(...
    'N',ones(1,cx)*rx,...%Observations
    'min',zeros(1,cx),...%Minimum values (zeroth quartile).
    'Q1',zeros(1,cx),... %First quartile.
    'med',zeros(1,cx),...%Median (second quartile).
    'Q3',zeros(1,cx),... %Third quartile.
    'max',zeros(1,cx),...%Maximum values (fourth quartile).
    'ran',zeros(1,cx),...%Range.
    'avg',zeros(1,cx),...%Arithmetic means.
    'std',zeros(1,cx),...%Standard deviation.
    'var',zeros(1,cx),...%Variance.
    'skw',zeros(1,cx),...%Skewness.
    'kur',zeros(1,cx));  %Kurtosis.
for ii = 1:cx
    X(:,ii) = sort(X(:,ii)); %CHEATING.
    
    %The average must be known first..
    for jj = 1:rx %Add-up all entries in the dataset.
        summary.avg(ii) = summary.avg(ii) + X(jj,ii);
    end
    if isnan(summary.avg(ii)) == 1
        error('NaN detected in data set.');
    end
    summary.avg(ii) = summary.avg(ii)/rx; %Final average.
    
    %Standardized moments all hinge on the average.
    for jj = 1:rx
        buffer = X(jj,ii) - summary.avg(ii); %(X_{i}- E[X])
        buffer = buffer*buffer; %(X_{i}- E[X])^2
        summary.var(ii) = summary.var(ii) + buffer;
        buffer = buffer*buffer; %(X_{i}- E[X])^3
        summary.skw(ii) = summary.skw(ii) + buffer;
        buffer = buffer*buffer; %(X_{i}- E[X])^4
        summary.kur(ii) = summary.kur(ii) + buffer;
    end
    summary.var(ii) = summary.var(ii)/rx; %Final variance.
    summary.skw(ii) = summary.skw(ii)/(rx*summary.var(ii)^(3/2)); %Final skewness.
    summary.kur(ii) = summary.kur(ii)/(rx*summary.var(ii)^2); %Final kurtosis.
    summary.std(ii) = sqrt(summary.var(ii));%Standard deviation.
    
    
    summary.min(ii) = X(1,ii);%0th-quartile.
    summary.max(ii) = X(rx,ii);%100th-quartile.
    summary.ran(ii) = X(rx,ii) - X(1,ii);%Range.
    %50th-percentile (median).
    if mod(rx,2) == 0 %Even number of observations.
        half = rx/2;
        summary.med(ii) = (X(half,ii) + X(half + 1,ii))/2;
    else %Odd number of observations.
        half = (rx+1)/2;
        summary.med(ii) = X(half,ii);
    end
    
    %25th and 75th percentiles (quartiles)
    if mod(half,2) == 0 %Even # of observations.
        summary.Q1(ii) = X(half/2,ii);
        summary.Q3(ii) = X(half + half/2,ii);
    else %Odd # of observations.
        summary.Q1(ii) = X((half+1)/2,ii);
        summary.Q3(ii) = X(round((rx-half+1)/2),ii); %%??????
    end
        
end
end