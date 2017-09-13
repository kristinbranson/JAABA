function [mang, stdang, mepe]=flowAngErr(tu, tv, u, v, bord, varargin)
% return the Barron et al angular error.  bord is the pixel width of the
% border to be ingnored.
%
%   Author:  Stefan Roth, Department of Computer Science, TU Darmstadt
%   Contact: sroth@cs.tu-darmstadt.de
%   $Date: 2007-03-27 14:09:11 -0400 (Tue, 27 Mar 2007) $
%   $Revision: 252 $

% Copyright 2004-2007, Brown University, Providence, RI. USA
% Copyright 2007-2010 TU Darmstadt, Darmstadt, Germany.
% 
%                          All Rights Reserved
% 
% All commercial use of this software, whether direct or indirect, is
% strictly prohibited including, without limitation, incorporation into in
% a commercial product, use in a commercial service, or production of other
% artifacts for commercial purposes.     
%
% Permission to use, copy, modify, and distribute this software and its
% documentation for research purposes is hereby granted without fee,
% provided that the above copyright notice appears in all copies and that
% both that copyright notice and this permission notice appear in
% supporting documentation, and that the name of the author and Brown
% University not be used in advertising or publicity pertaining to
% distribution of the software without specific, written prior permission.        
%
% For commercial uses contact the Technology Venture Office of Brown University
% 
% THE AUTHOR AND BROWN UNIVERSITY DISCLAIM ALL WARRANTIES WITH REGARD TO
% THIS SOFTWARE, INCLUDING ALL IMPLIED WARRANTIES OF MERCHANTABILITY AND
% FITNESS FOR ANY PARTICULAR PURPOSE.  IN NO EVENT SHALL THE AUTHOR OR
% BROWN UNIVERSITY BE LIABLE FOR ANY SPECIAL, INDIRECT OR CONSEQUENTIAL
% DAMAGES OR ANY DAMAGES WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR
% PROFITS, WHETHER IN AN ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS
% ACTION, ARISING OUT OF OR IN CONNECTION WITH THE USE OR PERFORMANCE OF
% THIS SOFTWARE.        

smallflow=0.0;

% if length(varargin) == 1
%     % mask
%     tu(~varargin{1}) = 0;
%     tv(~varargin{1}) = 0;
% end;
    
stu=tu(bord+1:end-bord,bord+1:end-bord);
stv=tv(bord+1:end-bord,bord+1:end-bord);
su=u(bord+1:end-bord,bord+1:end-bord);
sv=v(bord+1:end-bord,bord+1:end-bord);

% ignore a pixel if both u and v are zero 
%ind2=find((stu(:).*stv(:)|sv(:).*su(:))~=0);
ind2=find(abs(stu(:))>smallflow|abs(stv(:)>smallflow)); 
%length(ind2)
n=1.0./sqrt(su(ind2).^2+sv(ind2).^2+1);
un=su(ind2).*n;
vn=sv(ind2).*n;
tn=1./sqrt(stu(ind2).^2+stv(ind2).^2+1);
tun=stu(ind2).*tn;
tvn=stv(ind2).*tn;
ang=acos(un.*tun+vn.*tvn+(n.*tn));
mang=mean(ang);
mang=mang*180/pi;

if nargout >= 2
    stdang = std(ang*180/pi);
end;

if nargout == 3
    epe = sqrt((stu-su).^2 + (stv-sv).^2);
    epe = epe(ind2);
    mepe = mean(epe(:));
end;

% show which pixels were ignored
% tmp=zeros(size(su));
% tmp(ind2)=1;
% imagesc(tmp)