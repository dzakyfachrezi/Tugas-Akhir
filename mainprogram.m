%% MAIN PROGRAM FETAL ECG EXTRACTION WITH LMS FOR FETAL HEART RATE DETECTION
    % Designed & developed by Dzaky Fachrezi (18318040)
    % This program was used for final project in undergaduate program of
    % biomedical engineering, school of electric and informatics engineering,
    % bandung institute of technology

%% Dokumentasi
    % Input utama : 1 data sinyal dari dataset NIFEADB
    % Output utama : detak jantung janin (fHR) setiap 10 detik tersedia pada
    % variabel output globalresdat, data snr pra-ekstraksi (globalsnrpre), data
    % snr pasca-ekstraksi (globalsnrpost)
    % Keterangan kolom" pada globalresdat:
        % 1 : snr pasca ekstraksi untuk segmen residu kanal terpilih
        % 2 : jumlah fetal R anotasi per segmen
        % 3 : rata-rata interval R-R referensi
        % 4 : jumlah fetal R terdeteksi (sebelum validasi interval R-R)
        % 5 : rata-rata interval R-R terdeteksi (sebelum validasi interval R-R)
        % 6 : FP (sebelum validasi interval R-R)
        % 7 : FN (sebelum validasi interval R-R)
        % 8 : absolute error rata-rata interval R-R (sebelum validasi interval R-R)
        % 9 : jumlah deteksi missed beat
        % 10 : jumlah false detection
        % 11 : jumlah fetal R terdeteksi (setelah validasi interval R-R)
        % 12 : rata - rata interval R-R terdeteksi (setelah validasi interval R-R)
        % 13 : FP (setelah validasi interval R-R)
        % 14 : FN (setelah validasi interval R-R)
        % 15 : absolute error rata-rata interval R-R (setelah validasi interval R-R)
        % 16 : fHR referensi per 10 detik
        % 17 : fHR deteksi per 10 detik
    
    % no_seg : (jumlah segmen berdurasi 10 detik)

%% Read Signal & Input Parameters
% Sesuaikan parameter input dan IMPORT DATA ANOTASI FETAL-R DENGAN FORMAT
% MATRIX DAN NAMA: globalfQRS
dat = 'NR_10';  % sesuaikan dengan nama data sinyal (PAKAI UNDERSCORE CONTOH : 'ARR_01')
[signal,Fs,globaltm] = rdsamp('NR_10'); % sesuaikan dengan nama data sinyalnya (PAKAI UNDERSCORE CONTOH : 'ARR_01')
ch_annotate = 2; % kanal yang digunakan untuk anotasi (lihat data_info.xlsx)
f_ch = 1; % kanal terpilih (sesuai rata-rata SNR post ekstraksi tertinggi)
no_seg = 59; % jumlah segmen (lihat info_data.xlsx)

% Assign the maternal chest ECG (mECG) and abdominal ECG (aECG)
m_ch = 1;
a_ch = [2,3,4,5,6];
globalmECG = signal(:,m_ch);
globalaECG = signal(:,a_ch);

%% Signal Conditioning I (Baseline Wander & Powerline Interference Noise Removal)
% Zero-padding the signal
globalaECG_zp = [zeros(10000,length(a_ch));globalaECG;zeros(10000,length(a_ch))];
globalmECG_zp = [zeros(10000,length(m_ch));globalmECG;zeros(10000,length(m_ch))];

% Initialize variables
globalaECG_bwr = zeros(size(globalaECG_zp)); % abdominal ECG after baseline wander removal
globalaECG_pp = globalaECG_bwr;              % abdominal ECG after complete pre-processing
globalmECG_bwr = zeros(size(globalmECG_zp)); % maternal ECG after baseline wander removal
globalmECG_pp = globalmECG_bwr;              % maternal ECG after complete pre-processing

% Design HPF to remove baseline wander
Fstop1 = 0.5; % stopband frequency
Fpass1 = 1;   % passband frequency
Apass1 = 0.1; % passband ripple
Astop1 = 20;  % stopband attenuation

hpf = designfilt('highpassiir','StopbandFrequency',Fstop1,'PassbandFrequency',Fpass1,...
    'StopbandAttenuation',Astop1,'PassbandRipple',Apass1,'DesignMethod','butter',...
    'MatchExactly','stopband','SampleRate',Fs);
%fvtool(hpf)

% Design LPF to remove powerline interference
Fpass2 = 48;  % stopband frequency
Fstop2 = 50;  % passband frequency
Apass2 = 0.1; % passband ripple
Astop2 = 20;  % stopband attenuation

lpf = designfilt('lowpassiir','PassbandFrequency',Fpass2,'StopbandFrequency',Fstop2,...
    'PassbandRipple',Apass2,'StopbandAttenuation',Astop2,'DesignMethod','butter',...
    'MatchExactly','stopband','SampleRate',Fs);
%fvtool(lpf)

% Implement HPF & LPF with zero-phase digital filtering
for j = 1:length(a_ch)
    globalaECG_bwr(:,j) = filtfilt(hpf,globalaECG_zp(:,j));
    globalaECG_pp(:,j) = filtfilt(lpf,globalaECG_bwr(:,j));
end
for j = 1:length(m_ch)
    globalmECG_bwr(:,j) = filtfilt(hpf,globalmECG_zp(:,j));
    globalmECG_pp(:,j) = filtfilt(lpf,globalmECG_bwr(:,j));
end

% Trim the zero-padded signal
globalaECG_bwr = globalaECG_bwr(10001:length(globalaECG) + 10000,:);
globalmECG_bwr = globalmECG_bwr(10001:length(globalmECG) + 10000,:);
globalaECG_pp = globalaECG_pp(10001:length(globalaECG) + 10000,:);
globalmECG_pp = globalmECG_pp(10001:length(globalmECG) + 10000,:);

clear globalaECG_zp globalmECG_zp Fstop1 Fstop2 Fpass1 Fpass2 Astop1 Astop2 Apass1 Apass2 hpf lpf;

%% Signal Conditioning II (Alligning Maternal R-Peaks in Chest Lead and Abdominal Leads)
% Initialize variables
optimal_lag = zeros(size(a_ch)); % optimum lag of each signals to the mECG
aECG_alg = cell(size(a_ch)); % cell containing alligned aECG

for i = 1:length(a_ch)
    % find the time lag using maximum cross correlation
    [acor,lag] = xcorr(globalmECG_pp,globalaECG_pp(:,i));
    [~, I] = max(abs(acor));
    optimal_lag(i) = lag(I);

    % allign both signals based on optimal lag
    if optimal_lag(i) > 0
        addval = zeros(optimal_lag(i),1);
        aECG_alg{i} = [addval;globalaECG_pp(:,i)];
    else
        aECG_alg{i} = globalaECG_pp(-optimal_lag(i)+1:end,i);
    end
end

% Trim the edge difference
st_alg = max(optimal_lag);
if st_alg < 0
    st_alg = 0;
end

en_alg = -min(optimal_lag);
if en_alg < 0
    en_alg = 0;
end

len_alg = length(globalmECG_pp) - st_alg - en_alg;

globalaECG_pp = zeros(len_alg,length(a_ch));
for i = 1:length(a_ch)
    globalaECG_pp(:,i) = aECG_alg{i}(st_alg+1:st_alg+len_alg);
end

globalmECG_pp = globalmECG_pp(st_alg+1:st_alg+len_alg);

% fQRS reference correction
globalfQRS = globalfQRS + optimal_lag(ch_annotate) - st_alg;

clear acor lag I addval;

%% mECG Orientation Correction
% Find the polarity of mECG (1 : positive, 0 : negative)
y = sort(globalmECG_pp);
pol_mECG = mean(abs(y(end-floor(0.1*length(y))+1:end))) > mean(abs(y(1:floor(0.1*length(y)))));

% Orientation correction for mECG
if pol_mECG == 0
    globalmECG = -globalmECG;
    globalmECG_bwr = -globalmECG_bwr;
    globalmECG_pp = -globalmECG_pp;
end

clear y;

%% QRS Detection of mECG
% Maternal R-peaks detection with Pan-Tompkins Algorithm
[~,mqrs_i_raw,mdelay]=pan_tompkin(globalmECG_pp,Fs,0);
globalmR = zeros(length(mqrs_i_raw),2);
globalmR(:,1) = mqrs_i_raw; % maternal R-peaks index location
globalmR(:,2) = (1:length(mqrs_i_raw)); % indexing amount of mR

% Q and S detection of mECG
globalmQ = zeros(size(globalmR)); % maternal Q
globalmS = zeros(size(globalmR)); % maternal S
search_range = Fs*0.05; % search range of Q and S in sample unit (50 ms before and after R)

% Find the Q and S of each R-peaks
for i = 1:length(globalmR)
    if globalmR(i,1) <= search_range
        [~,tmp] = min(globalmECG_pp(1:globalmR(i,1)));
        globalmQ(i,1) = tmp;
        [~,tmp] = min(globalmECG_pp(globalmR(i,1):globalmR(i,1)+search_range));
        globalmS(i,1) = tmp + globalmR(i,1) - 1;
    elseif globalmR(i,1) + search_range >= length(globalmECG_pp)
        [~,tmp] = min(globalmECG_pp(globalmR(i,1)-search_range:globalmR(i,1)));
        globalmQ(i,1) = tmp + (globalmR(i,1) - search_range) - 1;
        [~,tmp] = min(globalmECG_pp(globalmR(i,1):end));
        globalmS(i,1) = tmp + globalmR(i,1) - 1;
    else
        [~,tmp] = min(globalmECG_pp(globalmR(i,1)-search_range:globalmR(i,1)));
        globalmQ(i,1) = tmp + (globalmR(i,1) - search_range) - 1;
        [~,tmp] = min(globalmECG_pp(globalmR(i,1):globalmR(i,1)+search_range));
        globalmS(i,1) = tmp + globalmR(i,1) - 1;
    end
end

globalmQ(:,2) = (1:length(globalmQ));
globalmS(:,2) = (1:length(globalmS));

% Check the validity of mQ and mS (mQ(i) < mS(i) and mS(i) < mQ(i+1)
% size of mQ, mR, and mS should be similar
valid = 0;
for i = 1:length(globalmR)
    if i == length(globalmR)
        if globalmQ(i,1) < globalmS(i,1)
            valid = 1;
        else
            valid = 0;
            break
        end
    else
        if globalmQ(i,1) < globalmS(i,1) && globalmS(i,1) < globalmQ(i+1,1)
            valid = 1;
        else
            valid = 0;
            break
        end
    end
end

if length(globalmQ) ~= length(globalmR) || length(globalmR) ~= length(globalmS)
    valid = 0;
    disp('the size of mQ, mR, and mS are unequal!')
end

if valid == 0
    disp('mQRS are not valid!')
end

clear mqrs_i_raw tmp valid search_range;

%% aECG Orientation Correction
pol_aECG = zeros(size(a_ch)); % aECG polarity matrix

for i = 1:length(a_ch)
    % find the polarity of aECG (1 : positive, 0 : negative)
    pospeaks = 0;
    negpeaks = 0;
    for j = 1:length(globalmR)
        if globalaECG_pp(globalmR(j,1),i) >= 0
            pospeaks = pospeaks + 1;
        else
            negpeaks = negpeaks + 1;
        end
    end
    
    if pospeaks > negpeaks
        pol_aECG(i) = 1;
    else
        pol_aECG(i) = 0;
    end

    % orientation correction for aECG
    if pol_aECG(i) == 0
        globalaECG(:,i) = -globalaECG(:,i);
        globalaECG_bwr(:,i) = -globalaECG_bwr(:,i);
        globalaECG_pp(:,i) = -globalaECG_pp(:,i);
    end
end

clear pospeaks negpeaks;

%% Segment for Testing (every 10s)
startPt = 0; % starting index - 1
seg_duration = Fs*10; % equals to 10 seconds
endPt = startPt + seg_duration; % ending index
globalsnrpre = zeros(no_seg,length(a_ch)); % snr pre-extraction-> column:channel ; row:segment
globalsnrpost = globalsnrpre; % snr post-extraction-> column:channel ; row:segment
globalart = globalsnrpre; % detection of potential artifacts in aECG segment all channels
globalresdat = zeros(no_seg,17); % summary of result
confmat = zeros(no_seg,6); % confusion matrix
globalfr_val = []; % fetal R index validated
globalrrdet = []; % RR interval validated
globalrrdetidx = [];
globalrrref = []; % RR interval reference
globalrrrefidx = [];
iterseg = 1; % number of iteration

while endPt <= length(globalmECG_pp)
    range_test = (startPt+1:endPt); % range index
    
    % Adjust the range of signals and time for each segments
    aECG = globalaECG(range_test,:);
    aECG_bwr = globalaECG_bwr(range_test,:);
    aECG_pp = globalaECG_pp(range_test,:);
    mECG = globalmECG(range_test);
    mECG_bwr = globalmECG_bwr(range_test);
    mECG_pp = globalmECG_pp(range_test);
    tm = globaltm(range_test);
    
    % Adjust the range of mQRS and fQRS
    mQ = globalmQ(globalmQ(:,1) >= (startPt+1) & globalmQ(:,1) <= endPt,:);
    mR = globalmR(globalmR(:,1) >= (startPt+1) & globalmR(:,1) <= endPt,:);
    mS = globalmS(globalmS(:,1) >= (startPt+1) & globalmS(:,1) <= endPt,:);
    mQ(:,1) = mQ(:,1) - startPt;
    mR(:,1) = mR(:,1) - startPt;
    mS(:,1) = mS(:,1) - startPt;
    fQRS = globalfQRS(globalfQRS >= (startPt+1) & globalfQRS <= endPt);
    fQRS = fQRS - startPt;
    
    %% Potential Artifacts Detection
    % Initialize variables
    art_aECG = zeros(size(a_ch)); % artifact flag (1 : detected, 0 : undetected)
    Psub = zeros(length(a_ch),5); % matrix of aECG sub_segment mean power
    
    % Measure the mean power for each aECG sub_segment (1 sub_seg = 1/5 length of aECG)
    for i = 1:length(a_ch)
        for j = 1:5
            if j == 5
                sub_seg = aECG_pp((j-1)*floor(length(aECG_pp(:,i))/5)+1:end,i);
            else
                sub_seg = aECG_pp((j-1)*floor(length(aECG_pp(:,i))/5)+1:j*floor(length(aECG_pp(:,i))/5),i);
            end
            sum_pow = sum(abs(sub_seg).^2);
            Psub(j,i) = sum_pow/length(sub_seg);
        end
    end
    
    % Check whether potential artifacts are found
    for i = 1:length(a_ch)
        for j = 1:5
            if Psub(j,i) > 1.5*(median(Psub(:,i)))
                art_aECG(i) = 1;
            end
        end
        if art_aECG(i) == 1
            %disp(['potential artifacts are found in aECG channel ',num2str(i)])
        end
    end
    
    globalart(iterseg,:) = art_aECG;

    clear sub_seg sum_pow Psub;
    % aECG with potential artifacts found are not included in the test

    %% mECG Component Decomposition from aECG with LMS
    % Initialize variables
    mECG_ab = zeros(size(aECG_pp)); % mECG component in abdominal channels estimated by the LMS
    mu = zeros(size(a_ch)); % step size
    nC = zeros(size(a_ch)); % adaptive filter length
    MSE = zeros(size(a_ch)); % MSE of each channel
    ext_row = Fs*5; % extra row for adapting the weight
    
    % Run the LMS
    for j = 1:length(a_ch)
        [mECG_ab(:,j),MSE(j),mu(j),nC(j)] = LMS_extraction([mECG_pp(ext_row:-1:1);mECG_pp],[aECG_pp(ext_row:-1:1,j);aECG_pp(:,j)],ext_row);
    end
    
    %% Compute Residual
    % Initialize variable
    residu = zeros(size(aECG_pp)); % extracted fECG from aECG
    
    % Subtract the mECG component from aECG and cancel out the mQRS
    for j = 1:length(a_ch)
        amp_aECG = zeros(size(mR)-2);
        amp_mECG_ab = amp_aECG;
        for i = 2:(length(mR) - 1)
            amp_aECG(i-1) = max(abs(aECG_pp(mR(i),j) - aECG_pp(mQ(i),j)),abs(aECG_pp(mR(i),j) - aECG_pp(mS(i),j)));
            amp_mECG_ab(i-1) = max(abs(mECG_ab(mR(i),j) - mECG_ab(mQ(i),j)),abs(mECG_ab(mR(i),j) - mECG_ab(mS(i),j)));
        end
        c1 = mean(amp_aECG);
        c2 = mean(amp_mECG_ab);
        % adjust the amplitude of mECG component
        mECG_ab(:,j) = mECG_ab(:,j)*(c1/c2);

        % substract to get residual
        residu(:,j) = aECG_pp(:,j) - mECG_ab(:,j);
    
        % cancel out mQRS
        elim_range = Fs*0.05;
        for i = 1:length(mR)
            if mR(i,1) <= elim_range
                residu(1:mR(i,1) + elim_range,j) = 0;
            elseif mR(i,1) + elim_range >= length(residu(:,j))
                residu(mR(i,1) - elim_range:end,j) = 0;
            else
                residu(mR(i,1) - elim_range:mR(i,1) + elim_range,j) = 0;
            end
        end
    end
    
    clear amp_mECG_ab amp_aECG c1 c2 elim_range;

    % correction of residual orientation (this process is just for
    % visualization necessity, the change of residuals orientation wont
    % disturb the fetal R peaks detection
    pol_res = adjustrespol(dat);
    for i = 1:length(a_ch)
        if pol_res(i) == 0
            residu(:,i) = -residu(:,i);
        end
    end
    
    %% SNR Measurement
    % Initialize variables
    sigpre = cell(size(a_ch));  % signal (QRS complex)
    noipre = cell(size(a_ch));  % noise (the rest)
    snrpre = zeros(size(a_ch)); % signal to noise ratio for each residual
    sigpost = cell(size(a_ch));  % signal (QRS complex)
    noipost = cell(size(a_ch));  % noise (the rest)
    snrpost = zeros(size(a_ch)); % signal to noise ratio for each residual
    
    % Seperate the signal (QRS complex) from the noise (the rest) and check the
    % validity (sum size of signal and noise should be equal to the size of
    % residual
    idx = (1:length(residu))';
    idx_sig = [];
    sig_range = Fs*0.04;

    for i = 1:length(fQRS)
        if fQRS(i) <= sig_range
            tmp = (1:fQRS(i)+sig_range);
        elseif fQRS(i) + sig_range >= length(residu)
            tmp = (fQRS(i)-sig_range:length(residu));
        else
            tmp = (fQRS(i)-sig_range:fQRS(i)+sig_range);
        end
        idx_sig = [idx_sig;tmp'];
    end
    
    idx_sig = unique(idx_sig);
    idx_sig = sort(idx_sig);
    elmt = ismember(idx,idx_sig);
    idx_noi = idx(~elmt);
    
    % Check the validity of signals and noises size
    cmb = sort([idx_sig;idx_noi]);
    valid = 1;
    for i = length(idx)
        if cmb(i) ~= idx(i)
            valid = 0;
            break
        end
    end
    
    if valid == 0
        disp('the sum size of sig and noise unequals to residu')
    end
    
    % Measure SNR_post extraction
    for i = 1:length(a_ch)
        sigpost{i} = residu(idx_sig,i);
        noipost{i} = residu(idx_noi,i);
        sig_pow = sum(abs(sigpost{i}).^2);
        noi_pow = sum(abs(noipost{i}).^2);
        Pf = sig_pow/length(sigpost{i});
        Pn = noi_pow/length(noipost{i});
        snrpost(i) = 10*log10(Pf/Pn);
    end
    
    % Measure SNR_pre extraction
    for i = 1:length(a_ch)
        sigpre{i} = aECG_pp(idx_sig,i);
        noipre{i} = aECG_pp(idx_noi,i);
        sig_pow = sum(abs(sigpre{i}).^2);
        noi_pow = sum(abs(noipre{i}).^2);
        Pf = sig_pow/length(sigpre{i});
        Pn = noi_pow/length(noipre{i});
        snrpre(i) = 10*log10(Pf/Pn);
    end

    % Assign to globalsnr and globalpol
    globalsnrpre(iterseg,:) = snrpre;
    globalsnrpost(iterseg,:) = snrpost;

    clear idx idx_sig sig_range idx_noi tmp elmt cmb valid sig_pow noi_pow Pf Pn;
    
    %% Residual Channel Selection
    % based on the highest average snrpost channels after all snr data is
    % collected
    fECG = residu(:,f_ch);
    
    %% R-Peaks Detection for fECG
    % Run Pan-Tompkins algorithm to get the R-Peaks of fECG (selected residual)
    [~,~,fqrs_i_raw,fdelay]=pan_tompkin_fecg(fECG,Fs,0);
    fR_det = fqrs_i_raw'; % Fetal R detected (unvalidated)
    
    clear fqrs_i_raw;
    
    % Confusion matrix of fetal r detection
    [TP_pre,FP_pre,FN_pre] = confmatrix_fR(fQRS,fR_det,length(fECG),Fs);

    %% R-R Interval Validation
    % Assign detected R-R interval
    rawRR = zeros(length(fR_det)-1,1);
    
    for i = 1:length(fR_det) - 1
        rawRR(i) = fR_det(i+1) - fR_det(i);
    end
    
    % Missed beat
    RR_mb = [];
    % check the first fR_det
    mbst = 0;
    RRst = fR_det(1) - 1;
    if RRst > 1.3*median(rawRR)
        intdivst = RRst / median(rawRR);
        if mod(RRst,median(rawRR)) < 0.6*median(rawRR)
            intdivst = floor(intdivst);
        else
            intdivst = ceil(intdivst);
        end
        tmpst = minQuotient(RRst,intdivst);
        if intdivst > 1 && sum(tmpst(1:end-1)) <= RRst
            RR_mb = [RR_mb;tmpst(1:end-1)];
            mbst = 1;
        end
    end
    % check the middle part of fR_det
    intdiv = [];
    for i = 1:length(rawRR)
        if rawRR(i) > 1.3*median(rawRR)
            intdiv(i) = rawRR(i) / median(rawRR);
            if mod(rawRR(i),median(rawRR)) < 0.6*median(rawRR)
                intdiv(i) = floor(intdiv(i));
            else
                intdiv(i) = ceil(intdiv(i));
            end
            tmp = minQuotient(rawRR(i),intdiv(i));
        else
            tmp = rawRR(i);
        end
        RR_mb = [RR_mb;tmp];
    end
    % check the last fR_det
    RRen = length(fECG) - fR_det(end);
    if RRen > 1.3*median(rawRR)
        intdiven = RRen / median(rawRR);
        if mod(RRen,median(rawRR)) < 0.6*median(rawRR)
            intdiven = floor(intdiven);
        else
            intdiven = ceil(intdiven);
        end
        tmpen = minQuotient(RRst,intdiven);
        if intdiven > 1 && sum(tmpen(1:end-1)) <= RRen
            RR_mb = [RR_mb;tmpen(1:end-1)];
        end
    end
    num_mb = length(RR_mb) - length(rawRR);
    
    % False detection
    RR_fd = [];
    for i = 1:length(RR_mb)
        if RR_mb(i) >= (Fs*0.2)
            RR_fd = [RR_fd;RR_mb(i)];
        end
    end
    num_fd = length(RR_mb) - length(RR_fd);
    
    % Misplaced
    RR_mp = RR_fd;
    for i = 1:length(RR_mp) - 1
        if (RR_mp(i)<0.9*median(RR_mp) && RR_mp(i+1)>1.1*median(RR_mp)) || (RR_mp(i+1)<0.9*median(RR_mp) && RR_mp(i)>1.1*median(RR_mp))
            tmp = minQuotient((RR_mp(i)+RR_mp(i+1)),2);
            RR_mp(i) = tmp(1);
            RR_mp(i+1) = tmp(2);
        end
    end
    
    % Validated R-R interval
    RRdet = RR_mp;
    
    % Final fetal R-peaks (validated)
    fR_val = zeros(length(RRdet)+1,1);
    % assign first fR_val
    if mbst == 0
        fR_val(1) = fR_det(1);
    else
        fR_val(1) = fR_det(1) - sum(tmpst(1:end-1));
    end
    % assign the rest of fR_val based on the RRdet
    for i = 1:length(RRdet)
        fR_val(i+1) = fR_val(i) + RRdet(i);
    end
    
    % Confusion matrix of fetal R detection after R-R validation
    [TP_post,FP_post,FN_post] = confmatrix_fR(fQRS,fR_val,length(fECG),Fs);

    clear RRst RRen intdivst intdiven tmpst tmpen intdiv tmp;
    
    %% fHR Calculation
    % Assign reference R-R interval
    RRref = zeros(length(fQRS)-1,1);
    for i = 1:length(fQRS)-1
        RRref(i) = fQRS(i+1) - fQRS(i);
    end
    
    % Calculate fHR
    mean_RRdet = mean(RRdet) * 1000 / Fs; % mean detected (validated) R-R interval
    mean_RRref = mean(RRref) * 1000 / Fs; % mean reference R-R interval
    mean_RRpre = mean(rawRR) * 1000 / Fs; % mean detected (unvalidated) R-R interval
    fHR_det = 60000 / mean_RRdet; % detected fHR (validated)
    fHR_ref = 60000 / mean_RRref; % reference fHR

    % Absolute & squared error
    AE_pre = abs(mean_RRref - mean_RRpre);
    AE_post = abs(mean_RRref - mean_RRdet);
    
    resdat = [snrpost(f_ch),length(fQRS),mean_RRref,length(fR_det),mean_RRpre,...
        FP_pre,FN_pre,AE_pre,num_mb,num_fd,length(fR_val),mean_RRdet,FP_post,...
        FN_post,AE_post,fHR_ref,fHR_det];

    tmp_confmat = [TP_pre,FP_pre,FN_pre,TP_post,FP_post,FN_post];

    globalresdat(iterseg,:) = resdat;
    confmat(iterseg,:) = tmp_confmat;
    globalfr_val = [globalfr_val;(fR_val + startPt)];
    globalrrdet = [globalrrdet;RRdet];
    globalrrdetidx = [globalrrdetidx;ones(length(RRdet),1)*iterseg];
    globalrrref = [globalrrref;RRref];
    globalrrrefidx = [globalrrrefidx;ones(length(RRref),1)*iterseg];
    
    startPt = startPt + (Fs*10);
    endPt = startPt + (Fs*10);
    iterseg = iterseg + 1;
end

% Untuk tabel 1
pretot_det_gagal = sum(confmat(:,2)) + sum(confmat(:,3));
preper_det_gagal = round(pretot_det_gagal*100/sum(globalresdat(:,2)),2);
pretot_mae = round(nanmean(globalresdat(:,8)),2);

postot_det_gagal = sum(confmat(:,5)) + sum(confmat(:,6));
posper_det_gagal = round(postot_det_gagal*100/sum(globalresdat(:,2)),2);
postot_mae = round(nanmean(globalresdat(:,15)),2);

tabel_1 = [sum(globalresdat(:,2)),pretot_det_gagal,preper_det_gagal,pretot_mae,...
    postot_det_gagal,posper_det_gagal,postot_mae];