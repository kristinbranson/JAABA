function y = getBoutLabel(jdata, expi, fly_id)
% Convert labels into the behavior bout label format used by the structured SVM code
    labels = jdata.GetLabels(expi,fly_id);
    [~, inds] = sort(labels.t0s);
    y = struct('num_bouts', [length(inds)], 'bouts', struct([]), 'fly_id', fly_id);
    for i=1:length(inds)
        bout = struct('start_frame', labels.t0s(inds(i)), ...
            'end_frame', labels.t1s(inds(i)),...  % TODO: assumes frame t1 not inside the bout (is this correct?)
            'behavior', labels.names{inds(i)});
        if i==1, y.bouts = bout; else y.bouts(i) = bout; end
    end
end

