myfile = 'soundrecord_985';
fc = 98.5e6;
frequencyTolerance = 10000;
packetTimeTolerance = 6.4e-4;
gapTime = 5e-3;
pkg load signal
graphics_toolkit('gnuplot')


for i = 5
  fs = 5e6;
  %system(sprintf('./signal_record.py %s %d %d', myfile, fs, fc)) %./signal_record.py must be in same directory as this run file.
  %system(sprintf('mkdir ~/Pictures/samp_%dMHz', i))
  for f = 10
    NFFT = 1024;
    [channels, thresh] = bandwidth_dur_estimator(myfile, fs, fc, NFFT, gapTime, packetTimeTolerance, frequencyTolerance);
    %saveas(1, sprintf('~/Pictures/samp_%dMHz/NFFT_%d_specgram.png', i, NFFT))
    %saveas(2, sprintf('~/Pictures/samp_%dMHz/NFFT_%d_histogram.png', i, NFFT))
  end
end