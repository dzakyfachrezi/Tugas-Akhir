function [TP,FP,FN] = confmatrix_fR(fQRS,fR_det,len_res)
    % Inisiasi variabel
    FN = 0;
    TP = 0;
    FP = 0;
    
    for i = 1:length(fQRS)
        % Rentang TP
        if fQRS(i) <= 62
            tmp = (1:fQRS(i)+62);
        elseif fQRS(i) + 62 >= len_res
            tmp = (fQRS(i)-62:len_res);
        else
            tmp = (fQRS(i)-62:fQRS(i)+62);
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