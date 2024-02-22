function pol_res = adjustrespol(dat)
    if strcmp(dat,'NR_02')
        pol_res = [1,1,0,0,0];
    elseif strcmp(dat,'NR_07')
        pol_res = [0,1,1,1,0];
    elseif strcmp(dat,'NR_10')
        pol_res = [0,0,0,0,0];
    elseif strcmp(dat,'NR_11')
        pol_res = [1,0,0,0,0];
    elseif strcmp(dat,'NR_14')
        pol_res = [0,1,1,1,1];
    elseif strcmp(dat,'ARR_01')
        pol_res = [1,1,0,0,0];
    elseif strcmp(dat,'ARR_03')
        pol_res = [0,1,0,0,0];
    elseif strcmp(dat,'ARR_09')
        pol_res = [1,0,1,1,1];
    elseif strcmp(dat,'ARR_11')
        pol_res = [1,0,0,0,0];
    elseif strcmp(dat,'ARR_12')
        pol_res = [0,0,0,0,0];
end