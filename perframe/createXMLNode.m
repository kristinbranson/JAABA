function node = createXMLNode(docNode,name,value)

node = docNode.createElement(name);

fs = fieldnames(value);
for ndx = 1:numel(fs)
  curVal = value.(fs{ndx});
  switch class(curVal)
    case 'struct'
      node.appendChild(createXMLNode(docNode,fs{ndx},value.(fs{ndx})));
    case 'cell'
      valStr = curVal{1};
      for vNdx = 2:numel(curVal)
        valStr = [valStr ',' curVal{vNdx}];
      end
      node.setAttribute(fs{ndx},valStr);
    case 'double'
      valStr = num2str(curVal(1));
      for vNdx = 2:numel(curVal)
        valStr = [valStr ',' num2str(curVal(vNdx))];
      end      
      node.setAttribute(fs{ndx},valStr);
    case 'single'
      valStr = num2str(curVal(1));
      for vNdx = 2:numel(curVal)
        valStr = [valStr ',' num2str(curVal(vNdx))];
      end      
      node.setAttribute(fs{ndx},valStr);
    case 'logical'
      valStr = num2str(curVal(1));
      for vNdx = 2:numel(curVal)
        valStr = [valStr ',' num2str(curVal(vNdx))];
      end      
      node.setAttribute(fs{ndx},valStr);
    case 'char'
      valStr = curVal;
      node.setAttribute(fs{ndx},valStr);
    otherwise
      fprintf('Unknown type');
  end
      
end