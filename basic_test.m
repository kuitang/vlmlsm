D = [2 5 6; 3 -2 4;  1 2 -7.2];
Vm(3,3,2) = 0;
Vm(:,:,1) = [1 2 4; 1 1 1; 4 1 0];
Vm(:,:,2) = [1 2 4; 2 3 4; 3 4 5];
Vi = [0 1 0; 1 0 2; 0 2 0];
W  = [0 1 0; 1 0 1; 0 1 0];

[xbf, ebf] = MultiLabelSubModularBruteForce(D, W, Vi, Vm)
[x, ebf]   = MultiLabelSubModular(D, W, Vi, Vm)
