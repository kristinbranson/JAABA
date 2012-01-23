function varargout = compute_closestfly_nose2ell_angle_30tomin30(trx,n,varargin)

anglerange = [30,-30]*pi/180;
varargout = compute_closestfly_nose2ell_anglerange(trx,n,anglerange,varargin{:});
