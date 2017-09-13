function url = lineReportURL(assayName,lineNames)
% assayName: one of {'box' 'bowl' ... }
% lineNames: a cellstr of line names.
%
% url: URL for informatics line report

assert(ischar(assayName),'assayName must be the name of an Olympiad assay.');
assert(iscellstr(lineNames),'lineNames must be a cell array of line names.');

URLBASE = 'http://informatics/cgi-bin/';

cgiStr = sprintf('olympiad_%s.cgi?Search=Search',assayName);

if ismac
    % On mac, the built-in web() shells out to unix directly and bash
    % thinks a raw semicolon is a command separator.
    lineStr = sprintf('\\;line=%s',lineNames{:});
else
    lineStr = sprintf(';line=%s',lineNames{:});
end

url = [URLBASE cgiStr lineStr];

end