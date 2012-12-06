function x = getBoutData(jdata, expi, fly_id)
    x = struct('fname', sprintf('%s/%s', jdata.expdirs{expi}, jdata.trxfilename), 'fly_id', fly_id);
end

