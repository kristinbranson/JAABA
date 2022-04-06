classdef file_object < handle
    % Simple file object that fcloses() the fid when no longer referenced.
    % Can be used just like a fid in most cases.
    
    properties (Access = protected)
        fid_ = []
        uid_ = []
        file_name_ = []
    end
    
    methods
        function self = file_object(varargin) 
            [fid, msg] = fopen(varargin{:}) ;
            if fid >= 0 ,
                self.fid_ = fid ;
                self.uid_ = randi([1 1000000]) ;
                self.file_name_ = varargin{1} ;
                %warning('In file_object with UID %d, on file %s, just fopened a file, got fid %d', self.uid_, self.file_name_, fid) ;
            else
                error('file_object:unable_to_open', msg) ;
            end
        end

        function result = fid(self)
            if ~isempty(self.fid_) && self.fid_>=0 ,
                result = self.fid_ ;
            else
                result = -1 ;
            end
        end
        
        function fprintf(self, varargin)            
            fprintf(self.fid(), varargin{:}) ;
        end
        
        function varargout = fread(self, varargin)
            if nargout==0 ,
                fread(self.fid(), varargin{:}) ;
            elseif nargout==1 ,
                A = fread(self.fid(), varargin{:}) ;
                varargout = {A} ;
            elseif nargout==2 ,
                [A,count] = fread(self.fid(), varargin{:}) ;
                varargout = {A, count} ;                
            else
                error('file_object:too_many_outputs', 'fread() can''t return more than two results') ;
            end
        end
        
        function result = ftell(self, varargin)
            result = ftell(self.fid(), varargin{:}) ;
        end
        
        function fseek(self, varargin)
            fseek(self.fid(), varargin{:}) ;
        end

        function result = lt(self, n)
            % Want fid<n comparisons to work properly
            result = (self.fid()<n) ;    
        end
        
        function result = lte(self, n)
            % Want fid<=n comparisons to work properly
            result = (self.fid()<=n) ;    
        end
        
        function result = gt(self, n)
            % Want fid>n comparisons to work properly
            result = (self.fid()>n) ;    
        end
        
        function result = gte(self, n)
            % Want fid>=n comparisons to work properly
            result = (self.fid()>=n) ;    
        end
        
        function fclose(self)
            fid = self.fid_ ;
            if isempty(fid) ,
                % make it the canonical empty array
                self.fid_ = [] ;
            else
                if fid >= 0 ,
%                     if isempty(self.uid_) ,
%                         warning('In file_object with *empty* UID, on file %s, about to fclose a file with fid %d', self.file_name_, fid) ;
%                     else
%                         warning('In file_object with UID %d, on file %s, about to fclose a file with fid %d', self.uid_, self.file_name_, fid) ;
%                     end
                    %try
                    fclose(fid) ;
%                     catch me
%                         if strcmp(me.identifier, 'MATLAB:FileIO:InvalidFid') ,
%                             % ignore it
%                         else
%                             rethrow(me) ;
%                         end
%                     end
                end
                % make it the canonical empty array
                self.fid_ = [] ;
            end
        end
        
        function delete(self)
            self.fclose() ;
        end
    end
end
