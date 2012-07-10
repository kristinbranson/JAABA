function SaveXMLParams(topNode,filename,name)

if nargin < 3,
  name = 'params';
end

docNode = com.mathworks.xml.XMLUtils.createDocument(name);
toc = docNode.getDocumentElement;
att = fieldnames(topNode);
for ndx = 1:numel(att)
  toc.appendChild(createXMLNode(docNode,att{ndx},topNode.(att{ndx})));
end
xmlwrite(filename,docNode);
