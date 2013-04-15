function varargout=unify()

controller=GrandlyUnifyController();

if nargout>=1 ,
  varargout{1}=controller;
end

end
