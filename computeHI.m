fHR = zeros(60,2); % kolom 1 untuk fHR referensi dan kolom 2 untuk fHR deteksi
% run terlebih dahulu hingga baris ini lalu input nilai-nilai fHR-nya

Tdecref = 0;
Tdecdet = 0;

for i=2:length(fHR)
    if fHR(i,1) < fHR(i-1,1)
        Tdecref = Tdecref + (1/6);
    end
end

for i=2:length(fHR)
    if fHR(i,2) < fHR(i-1,2)
        Tdecdet = Tdecdet + (1/6);
    end
end

HIref = round(Tdecref*100/min(fHR(:,1)),2);
HIdet = round(Tdecdet*100/min(fHR(:,2)),2);

res = [HIref,HIdet];