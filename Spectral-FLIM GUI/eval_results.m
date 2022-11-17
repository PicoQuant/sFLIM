clear all

name = 'test_results.mat';

load([name]);

amp = results.amp;
amp = reshape(amp,[size(amp,1)*size(tmp,2) size(amp,3)]);

m = max(max(amp));
m = 100*ceil(m/100);
b = 10:5:m;
histo_1 = [];

for a = 1:size(amp,3)
    histo(:,a) = mHist(amp(:,a),b);
    m1(a) = mean(amp(:,a));
    s1(a) = std(amp(:,a));
end

save([name(1:end-3) 'dat'],'histo','-ASCII')