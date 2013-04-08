function f = makeWrong(params)
% A "potentially wrong" problem just has a high L1 distance in marginals. 

    function tf = potentialWrong(trueNegBethe, t1M, approxLogZ, a1M)
        % An instance is wrong iff the one marginals are far (absolute
        % difference > params.margThresh) and the Bethe free energies of the
        % approximate and true solutions are also far (absolute differences
        % > params.betheThresh).        

        tf = any(abs(t1M - a1M)) > params.margThresh;
    end

    f = @potentialWrong;
end
