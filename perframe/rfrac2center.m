function [x1,y1,x2,y2] = rfrac2center(trx,fly,rfrac)

x1 = trx(fly).x_mm(1:end-1) + rfrac(1,:).*trx(fly).a_mm(1:end-1).*2.*cos(trx(fly).theta_mm(1:end-1)) - rfrac(2,:).*trx(fly).b_mm(1:end-1).*2.*sin(trx(fly).theta_mm(1:end-1));
y1 = trx(fly).y_mm(1:end-1) + rfrac(1,:).*trx(fly).a_mm(1:end-1).*2.*sin(trx(fly).theta_mm(1:end-1)) + rfrac(2,:).*trx(fly).b_mm(1:end-1).*2.*cos(trx(fly).theta_mm(1:end-1));
x2 = trx(fly).x_mm(2:end) + rfrac(1,:).*trx(fly).a_mm(2:end).*2.*cos(trx(fly).theta_mm(2:end)) - rfrac(2,:).*trx(fly).b_mm(2:end).*2.*sin(trx(fly).theta_mm(2:end));
y2 = trx(fly).y_mm(2:end) + rfrac(1,:).*trx(fly).a_mm(2:end).*2.*sin(trx(fly).theta_mm(2:end)) + rfrac(2,:).*trx(fly).b_mm(2:end).*2.*cos(trx(fly).theta_mm(2:end));
