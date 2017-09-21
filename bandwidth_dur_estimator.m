function [channels, thresh] = bandwidth_dur_estimator(myfile, fs, fc, NFFT, gapTolerance, packetTimeTolerance, frequencyTolerance)
%function bandwidth_dur_estimator reads in a complex binary file
close('all')
B = fs;
%NFFT = 1024;
FFTstep = NFFT/2;
y = read_complex_binary(myfile); %(using source --> Head --> file sink. files have no extension)

%fft processing
[S, c] = stft(y,NFFT,FFTstep,NFFT);
S = fftshift(S);

%creating frequency-axis vector for accurate plotting and indexing.
df = B/(2*NFFT);               
fOrig = 1:rows(S);
fNew = fOrig.*df + fc - B/2;  %centers fc to middle of graph and scales properly

%creating time-axis vector for accurate plotting and indexing.
dt = FFTstep / fs;
tOrig = 1:columns(S);
tNew = tOrig.*dt;          

%Display fft
LogS = 20*log10(abs(S));   %Decibel version of fft output.
figure(1); %subplot(2,1,1);
imagesc(tNew,fNew,LogS); axis('xy'); %colorbar;

Log = max(LogS,[],2);

[FF, XX] = hist(Log(:),250);
figure(2); hist(Log(:),250)

%Sgolay LPF
FF = sgolayfilt(FF,1);

%thresholding processing
%Method: find 
[m, ind] = findpeaks(1.01*max(FF) - FF); %find minima
m = FF(ind);
[mval maxind] = max(FF);
m(ind < maxind) = [];
ind(ind < maxind) = [];
if isempty(m)
  m = 0;
  ind = columns(FF);
  occurrenceThresh = ((mval - m(1)) * 0.3) + m(1);
else
  occurrenceThresh = ((mval - m(1)) * 0) + m(1);
end

X = maxind:ind(1);
Y = FF(X);
Xind = find((Y-occurrenceThresh) <= 0, 1);
thresh = XX(X(Xind)); %histogram x val between first peak and next min where hist y-val = occurencethresh


threshLog = LogS >= thresh; 
figure(3); %subplot(2,1,2); 
imagesc(tNew,fNew,threshLog); axis('xy')


%find start- and stop-frequencies of channels above threshold.
rFreqs = diff([1; ~any(threshLog,2); 1]);
startInd = find(rFreqs < 0);
endInd = find(rFreqs > 0)-1;
duration = endInd-startInd+1;
strInd = duration >= frequencyTolerance./df;
startInd = startInd(strInd);
endInd = endInd(strInd);

bandInd = [startInd'; endInd'];
if columns(bandInd) == 1
  bandvec = fNew(bandInd)';
else
  bandvec = fNew(bandInd);
end

%find start- and end-times of each packet. Use threshold values to filter noise.
channelCell = {};
[r,c] = size(bandInd);
bandvecLog = zeros(1,length(bandvec));
for i = 1:c
  channel = threshLog(bandInd(1,i):bandInd(2,i),:); 
  timeBand = any(channel,1);
  diffTime = diff([1, ~timeBand, 1]);
  startTimeInd = find(diffTime < 0);
  endTimeInd = find(diffTime > 0)-1;
  
  gapDurations = startTimeInd(2:end) - endTimeInd(1:end-1) + 1;
  gapStrInd = gapDurations <= gapTolerance/dt;
  startTimeInd([false, gapStrInd]) = [];
  endTimeInd([gapStrInd, false]) = [];
  
  timeDuration = endTimeInd - startTimeInd + 1;
  timeStrInd = timeDuration >= packetTimeTolerance/dt;
  startTimeInd = startTimeInd(timeStrInd);
  endTimeInd = endTimeInd(timeStrInd);
    
  if ~isempty(startTimeInd) && ~isempty(endTimeInd)
    channelCell{end+1} =  tNew([startTimeInd',endTimeInd']);
    bandvecLog(i) = true;
  end
end

bandvec = bandvec(:,logical(bandvecLog));
bandcell = {};

for i = bandvec
  bandcell{end+1} = i;
end

figure(1); hold on;
channels = struct('Bands',bandcell,'Packets',channelCell);
for i = 1:length(channels)bandvec = fNew(bandInd)';
  f = channels(i).Bands;
  tvec = channels(i).Packets;
  for j = 1:rows(tvec)
    t = tvec(j,:);
    plot(t([1 2 2 1 1]), f([1 1 2 2 1]),'k');
  end

end

%printf('found channels:\n')
%disp([channels.Bands])
endfunction