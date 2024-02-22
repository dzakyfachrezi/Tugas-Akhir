function [TP,FP,FN] = confmatrix_fR(fQRS,fR_det,len_res,Fs)
    if Fs == 1000
        search_range = 100;
    else
        search_range = 62;
    end
    % Inisiasi variabel
    FN = 0;
    TP = 0;
    FP = 0;
    
    for i = 1:length(fQRS)
        % Rentang TP
        if fQRS(i) <= search_range
            tmp = (1:fQRS(i)+search_range);
        elseif fQRS(i) + search_range >= len_res
            tmp = (fQRS(i)-search_range:len_res);
        else
            tmp = (fQRS(i)-search_range:fQRS(i)+search_range);
        end
    
        % FN & TP
        if sum(ismember(fR_det,tmp)) == 0
            FN = FN + 1;
        else
            TP = TP + 1;
        end
    end
    
    % FP
    FP = length(fR_det) - TP;
end