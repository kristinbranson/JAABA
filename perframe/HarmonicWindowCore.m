function res = HarmonicWindowCore(x,w,num_harmonic_curr)

r = (w-1)/2;
fil = cos(linspace(0,pi*num_harmonic_curr,w))/w*(num_harmonic_curr+1);
% res is nan from 1:r, N-r+1:N
res = myconv(x,fil,nan,'corr','same');
% boundary conditions: smallw used for frame
for smallr = 1:r-1,
  smallw = min(2*smallr+1,numel(x));
  fil = cos(linspace(0,pi*num_harmonic_curr,smallw))/smallw*(num_harmonic_curr+1);
  % start
  res(smallr+1) = sum(fil.*x(1:smallw));
  % end
  res(end-smallr) = sum(fil.*x(end-smallw+1:end));
end

