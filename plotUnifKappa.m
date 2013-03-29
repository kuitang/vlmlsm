function [ h, Ts, Ws, Es, Js, kappas ] = plotUnifKappa(nNodes, epsilon, ta, tb, wb, nPts)
% plotUnifKappa - Plot the kappa parameter for fully connected Ising
%
% h      - figure handle
% ta, tb - (uniform) theta limits
% wb     - (uniform) w upper limit
% nPts   - number of points along each side (total Npts^2 mesh points)

    thetas = linspace(ta, tb, nPts);
    ws     = linspace(epsilon, wb, nPts);    
    kappas = zeros(nPts, nPts);
    
    etas(nPts) = 0;
    js(nPts)   = 0;
    
    Wones  = ones(nNodes, nNodes);
    Wones  = Wones - diag(diag(Wones));
    tones  = ones(nNodes, 1);
    
    [Ts, Ws] = meshgrid(thetas, ws);
    
    for it = 1:nPts
        for iw = 1:nPts
            W     = sparse(ws(iw) * Wones);
            % Remove bias
            %
            % YOU CAN'T JUST KEEP REMOVING STUFF HERE! You really need to
            % meshgrid it in the beginning and change the 
            theta = thetas(it) * tones - 0.5 * sum(W, 2);
            Ts(it) = theta(1);            
            
            [j, e] = mooijParam(W, theta);
            % Should be uniform anyways
            js(iw)   = j(1,2);
            etas(it) = e(1);            
                        
            [A, B, ~] = BBP(theta, W, 1e-8, 1000);
            [intervalSz, kappas(it,iw), theoryBound] = getIntervalSz(A, B, W, epsilon);
        end
    end
    
    h = figure;
    contourf(Ts, Ws, log(kappas));
    colorbar;
    title('log(\kappa) vs. \theta and W');

    h = figure;
    [Es, Js] = meshgrid(etas, js);
    contourf(Es, Js, log(kappas));
    colorbar;
    title('log(\kappa) vs. \eta and J');    
end
