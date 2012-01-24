function [rfrac,isonfly] = center_of_rotation2(trx,fly,debug)
if ~exist('debug','var')
  debug = false;
end
N = 100;

% note that if dtheta = 0, this will produce 0

cost = cos(trx(fly).theta_mm);
sint = sin(trx(fly).theta_mm);
dacost = 2*diff(trx(fly).a_mm.*cost);
dbcost = 2*diff(trx(fly).b_mm.*cost);
dasint = 2*diff(trx(fly).a_mm.*sint);
dbsint = 2*diff(trx(fly).b_mm.*sint);
Z = dacost .* dbcost + dbsint.*dasint;
Minv = zeros(2,2,trx(fly).nframes-1);
Minv(1,1,:) = dbcost ./ Z;
Minv(1,2,:) = dbsint ./ Z;
Minv(2,1,:) = -dasint ./ Z;
Minv(2,2,:) = dacost ./ Z;
Minv = reshape(Minv,[4,trx(fly).nframes-1]);
x = [ diff(trx(fly).x_mm) ; diff(trx(fly).y_mm) ];
rfrac = -[Minv(1,:) .* x(1,:) + Minv(3,:) .* x(2,:);...
  Minv(2,:) .* x(1,:) + Minv(4,:) .* x(2,:)];
% when no rotation, set center of rotation to middle of fly by default
rfrac(isnan(rfrac)) = 0;
rfrac0 = rfrac;
isoutofbounds = sum(rfrac.^2,1) > 1;
psi = linspace(0,2*pi,N)';
cospsi = cos(psi);
sinpsi = sin(psi);
idx = find(isoutofbounds);
nout = length(idx);
if nout > 0,
  x1 = repmat(trx(fly).x_mm(idx),[N,1]) + repmat(2*trx(fly).a_mm(idx).*cost(idx),[N,1]).*repmat(cospsi,[1,nout]) - ...
    repmat(2*trx(fly).b_mm(idx).*sint(idx),[N,1]).*repmat(sinpsi,[1,nout]);
  y1 = repmat(trx(fly).y_mm(idx),[N,1]) + repmat(2*trx(fly).a_mm(idx).*sint(idx),[N,1]).*repmat(cospsi,[1,nout]) + ...
    repmat(2*trx(fly).b_mm(idx).*cost(idx),[N,1]).*repmat(sinpsi,[1,nout]);
  x2 = repmat(trx(fly).x_mm(idx+1),[N,1]) + repmat(2*trx(fly).a_mm(idx+1).*cost(idx+1),[N,1]).*repmat(cospsi,[1,nout]) - ...
    repmat(2*trx(fly).b_mm(idx+1).*sint(idx+1),[N,1]).*repmat(sinpsi,[1,nout]);
  y2 = repmat(trx(fly).y_mm(idx+1),[N,1]) + repmat(2*trx(fly).a_mm(idx+1).*sint(idx+1),[N,1]).*repmat(cospsi,[1,nout]) + ...
    repmat(2*trx(fly).b_mm(idx+1).*cost(idx+1),[N,1]).*repmat(sinpsi,[1,nout]);
  d = (x1 - x2).^2 + (y1 - y2).^2;
  [~,j] = min(d,[],1);
  rfrac(1,idx) = cospsi(j);
  rfrac(2,idx) = sinpsi(j);
end
isonfly = ~isoutofbounds;

if debug,
  ntry = 50;
  psitry = linspace(0,2*pi,ntry+2);
  psitry = psitry(2:end-1);
  [rhotry,psitry] = meshgrid(linspace(0,1,ntry),psitry);
  rtry = [rhotry(:).*cos(psitry(:)),rhotry(:).*sin(psitry(:))]';
  rtry = [rtry,[0;0]];
  colorsplot = 'rg';
  colorsplot2 = 'mc';
  [~,order] = sort(-abs(modrange(diff(trx(fly).theta_mm),-pi,pi)));
  for i = order,
    u = [0,0]; v = [0,0];
    uout = [0,0]; vout = [0,0];
    utry = zeros(2,length(rtry)); vtry = zeros(2,length(rtry));
    clf; hold on;
    for j = [i,i+1],
      u(j-i+1) = trx(fly).x_mm(j) + rfrac(1,i)*trx(fly).a_mm(j)*2*cos(trx(fly).theta_mm(j)) - rfrac(2,i)*trx(fly).b_mm(j)*2*sin(trx(fly).theta_mm(j));
      v(j-i+1) = trx(fly).y_mm(j) + rfrac(1,i)*trx(fly).a_mm(j)*2*sin(trx(fly).theta_mm(j)) + rfrac(2,i)*trx(fly).b_mm(j)*2*cos(trx(fly).theta_mm(j));
      uout(j-i+1) = trx(fly).x_mm(j) + rfrac0(1,i)*trx(fly).a_mm(j)*2*cos(trx(fly).theta_mm(j)) - rfrac0(2,i)*trx(fly).b_mm(j)*2*sin(trx(fly).theta_mm(j));
      vout(j-i+1) = trx(fly).y_mm(j) + rfrac0(1,i)*trx(fly).a_mm(j)*2*sin(trx(fly).theta_mm(j)) + rfrac0(2,i)*trx(fly).b_mm(j)*2*cos(trx(fly).theta_mm(j));
      utry(j-i+1,:) = trx(fly).x_mm(j) + rtry(1,:)*trx(fly).a_mm(j)*2*cos(trx(fly).theta_mm(j)) - rtry(2,:)*trx(fly).b_mm(j)*2*sin(trx(fly).theta_mm(j));
      vtry(j-i+1,:) = trx(fly).y_mm(j) + rtry(1,:)*trx(fly).a_mm(j)*2*sin(trx(fly).theta_mm(j)) + rtry(2,:)*trx(fly).b_mm(j)*2*cos(trx(fly).theta_mm(j));
      ellipsedraw(trx(fly).a_mm(j)*2,trx(fly).b_mm(j)*2,trx(fly).x_mm(j),trx(fly).y_mm(j),trx(fly).theta_mm(j),colorsplot(j-i+1));
      plot(trx(fly).x_mm(j)+trx(fly).a_mm(j)*2*cos(trx(fly).theta_mm(j)),trx(fly).y_mm(j)+trx(fly).a_mm(j)*2*sin(trx(fly).theta_mm(j)),'x','color',colorsplot(j-i+1));
      htest(j-i+1) = plot([trx(fly).x_mm(j),u(j-i+1),uout(j-i+1)],[trx(fly).y_mm(j),v(j-i+1),vout(j-i+1)],'o-','color',colorsplot(j-i+1),'markerfacecolor',colorsplot(j-i+1)); %#ok<AGROW>
    end
    d = (utry(1,:)-utry(2,:)).^2 + (vtry(1,:)-vtry(2,:)).^2;
    [dexp,k] = min(d);
    uexp = utry(:,k);
    vexp = vtry(:,k);
    dexp = sqrt(dexp);
    dtest = sqrt(diff(u).^2 + diff(v).^2);
    dout = sqrt(diff(uout).^2 + diff(vout).^2);
    plot(uexp,vexp,'k.-');
    plot(u,v,'b.-');
    plot(uout,vout,'b.-');
    for j = [i,i+1],
      hexp(j-i+1) = plot([uexp(j-i+1),trx(fly).x_mm(j)],[vexp(j-i+1),trx(fly).y_mm(j)],'o-','color',colorsplot2(j-i+1),'markerfacecolor',colorsplot2(j-i+1)); %#ok<AGROW>
    end
    ax = nan(1,4);
    ax1 = nan(1,4);
    for j = [i,i+1],
      [ax1(1),ax1(2),ax1(3),ax1(4)] = ellipse_to_bounding_box(trx(fly).x_mm(j),trx(fly).y_mm(j),2*trx(fly).a_mm(j),2*trx(fly).b_mm(j),trx(fly).theta_mm(j));
      ax([1,3]) = min(ax([1,3]),ax1([1,3]));
      ax([2,4]) = max(ax([2,4]),ax1([2,4]));
    end
    ax = ax + [-3,3,-3,3];
    legend([htest,hexp],'analytic1','analytic2','empirical1','empirical2');
    title(sprintf('red = frame %d, green = frame %d, dtest = %f, dout = %f, dexp = %f',i,i+1,dtest,dout,dexp));
    axis equal
    axis(ax)
    input('');
  end
end
