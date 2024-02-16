function necg = normecg(ecg,start,stop)
% Normalize signal based on min-max scaling of the signal segment
%if size(ecg,2)>size(ecg,1)
   % ecg = ecg';
%end

try
    % min-max scaling
    scalings=1./(max(ecg(start:stop,:))-min(ecg(start:stop,:)));
    % detrend the signal
    rescaled_ecg = bsxfun(@times,ecg,scalings);
    detrendedInput = bsxfun(@minus,rescaled_ecg,mean(rescaled_ecg(start:stop,:)));
    % take the tanh of the signal to avoid outliers
necg = tanh(detrendedInput);
catch ME
    rethrow(ME);
    for enb=1:length(ME.stack); disp(ME.stack(enb)); end;
    necg = ecg;
end