function f = makeWrong(params)
    % Make a "wrong" function handle that parameterizes the thresholds

    function tf = wrong(trueNegBethe, t1M, approxLogZ, a1M)
        % An instance is wrong iff the one marginals are far (absolute
        % difference > params.margThresh) and the Bethe free energies of the
        % approximate and true solutions are also far (absolute differences
        % > params.betheThresh).        
        tf = false;

        % Recall that for approximate solutions, Bethe = -logZ.        
        if trueNegBethe - approxLogZ > params.betheThresh && ...
            any(abs(t1M - a1M)) > params.margThresh
                tf = true;
                return;
        end        
    end

    f = @wrong;
end
