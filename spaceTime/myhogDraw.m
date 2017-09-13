function V = myhogDraw( H, w , directed)

if nargin<3
  directed = false;
end

nFold=1; s=size(H); s(3)=s(3)/nFold; w0=H; H=zeros(s);
for o=0:nFold-1, H=H+w0(:,:,(1:s(3))+o*s(3)); end;

% construct a "glyph" for each orientaion
if(nargin<2 || isempty(w)), w=15; end
bar=zeros(w,w); 
if directed,
  bar(round(.45*w):round(.55*w),(round(0.6*w)):end)=1;
else
  bar(round(.45*w):round(.55*w),:)=1;
end
bars=zeros([size(bar) s(3)]);
for o=1:s(3), 
  if directed,
    bars(:,:,o)=imrotate(bar,-(o-1)*360/s(3),'crop'); 
  else
    bars(:,:,o)=imrotate(bar,-(o-1)*180/s(3),'crop'); 
  end
end

% make pictures of positive weights by adding up weighted glyphs
H(H<0)=0; V=zeros(w*s(1:2));
for r=1:s(1), rs=(1:w)+(r-1)*w;
  for c=1:s(2), cs=(1:w)+(c-1)*w;
    for o=1:s(3), V(rs,cs)=V(rs,cs)+bars(:,:,o)*H(r,c,o); end
  end
end
