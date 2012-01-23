function varargout = compute_closestfly_nose2ell_angle_min20to20(trx,n,varargin)

anglerange = [-20,20]*pi/180;
varargout = compute_closestfly_nose2ell_anglerange(trx,n,anglerange,varargin{:});
