function fileNameRel=fileNameRelFromAbs(fileNameAbs)

[~,name,ext]=fileparts(fileNameAbs);
fileNameRel=[name ext];

end
