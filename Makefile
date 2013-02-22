BetheApprox_mex.mexmaci64: BetheApprox_mex.cpp BetheApprox.cpp BetheApprox.h MinSum.cpp MinSum.h
	mex -largeArrayDims -g BetheApprox_mex.cpp MinSum.cpp BetheApprox.cpp graph.cpp maxflow.cpp

