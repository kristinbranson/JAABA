function gfdata = getframe_initialize(haxes)

gfdata.haxes = haxes;
hfig = get(haxes,'parent');
gfdata.hfig = hfig;
gfdata.units = get(hfig,'Units');
gfdata.pos = get(hfig,'Position');
pixelsperinch = get(0,'screenpixelsperInch');
pos =  get(hfig,'position');
% If units are normalized, then re-scale paperposition calculation
% to use pixel values based on the screen size.
if strcmp(get(hfig, 'units'), 'normalized')
  screensize = get(0, 'screensize');
  screenwidth = screensize(3);
  screenheight = screensize(4);
  pos = [(pos(1) * screenwidth + 1) (pos(2) * screenheight + 1) (pos(3) * screenwidth) (pos(4) * screenheight)];
end
set(hfig, 'paperposition', pos./pixelsperinch);
renderer = get(hfig,'renderer');
if strcmp(renderer,'painters')
  renderer = 'zbuffer';
end
%Turn off warning in case opengl is not supported and
%hardcopy needs to use zbuffer
warnstate = warning('off');
gfdata.warnstate = warnstate;
noanimate('save',get(haxes,'parent'));
gfdata.hardcopy_args = {haxes,['-d',renderer],['-r',num2str(round(pixelsperinch))]};

% TODO: put this in cleanup
% noanimate('restore',get(haxes,'parent'));
% warning(warnstate);