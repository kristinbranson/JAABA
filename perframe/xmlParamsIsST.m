function tf = xmlParamsIsST(filename)

DOMnode = xmlread(filename);
n = DOMnode.getDocumentElement();
c = n.getChildNodes;
len = c.getLength;
for i = 0:len-1
  if strcmp(c.item(i).getNodeName,'st')
    tf = true;
    return;
  end
end

tf = false;
  