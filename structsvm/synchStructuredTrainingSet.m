function ids=synchStructuredTrainingSet( jdata, fname )
    ex = struct('x',{},'y',{});
    m = 0;
    i = 1;
    for expi = 1:jdata.nexps, m = max(m,jdata.nflies_per_exp(expi)); end
    ids = zeros(jdata.nexps, m);
    
    for expi = 1:jdata.nexps,
        for fly = 1:jdata.nflies_per_exp(expi)
            ex(i) = struct('x', getBoutData(jdata, expi, fly), ...
                'y', getBoutLabel(jdata, expi, fly));
            
            query = struct('method', 'add_example', 'x', ex(i).x, 'y', ex(i).y);
        
            response = jsonrpc_request(jdata.h.ip_address, jdata.h.port, query);
            if isfield(response, 'error'),  
                warning(response.error); 
                return; 
            end
            if ~isfield(response, 'index'),
                warning('No "index" field found');
                return
            end
            
            ids(expi,fly) = response.index;
            i = i+1;
        end
    end
    
    if nargin > 1,
        fid = fopen(fname, 'w');
        fprintf(fout, '%s', savejson(ex));
        fclose(fid);
    end
end

