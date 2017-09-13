function [success ips] = genInterestPoints(moviename,varargin)

[ipnames,doaskleftright] = myparse(varargin,'ipnames',{'food','mouth','perch'},...
  'doaskleftright',true);

annObj = HandleObj; % handle obj to store annotation info
hFig = playfmf('moviename',moviename,'annDataObj',annObj,'ipnames',ipnames,'doaskleftright',doaskleftright);
waitfor(hFig);
if isempty(annObj.data)
  success = false;
  ips = [];
  return;
end

fprintf(1,'Interest points generated:\n');
disp(annObj.data);

ips = struct();
for i = 1:numel(ipnames),
  ips.(ipnames{i}) = annObj.data.(ipnames{i});
end
if doaskleftright
  ips.face = annObj.data.lr;
end
success = true;