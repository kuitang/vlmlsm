all: BetheApprox_opt_mex.mexmaci64 BetheApprox_debug_mex.mexmaci64 BBP_MK_opt_mex.mexmaci64 makePotential_mex.mexmaci64 MultiLabelSubModular_mex.mexmaci64

MultiLabelSubModular_mex.mexmaci64: MultiLabelSubModular_mex.cpp
	mex -largeArrayDims -O MultiLabelSubModular_mex.cpp MinSum.cpp BetheApprox.cpp graph.cpp maxflow.cpp -o MultiLabelSubModular_mex.mexmaci64

makePotential_mex.mexmaci64: makePotential_mex.cpp
	mex -largeArrayDims -O makePotential_mex.cpp -o makePotential_mex.mexmaci64

BBP_MK_opt_mex.mexmaci64: BBP_MK_mex.cpp BetheApprox.h BetheApprox.cpp
	mex -largeArrayDims -O BBP_MK_mex.cpp MinSum.cpp BetheApprox.cpp graph.cpp maxflow.cpp -o BBP_MK_opt_mex.mexmaci64

BetheApprox_opt_mex.mexmaci64: BetheApprox_mex.cpp BetheApprox.cpp BetheApprox.h MinSum.cpp MinSum.h
	mex -largeArrayDims -O BetheApprox_mex.cpp MinSum.cpp BetheApprox.cpp graph.cpp maxflow.cpp -o BetheApprox_opt_mex.mexmaci64

BetheApprox_debug_mex.mexmaci64: BetheApprox_mex.cpp BetheApprox.cpp BetheApprox.h MinSum.cpp MinSum.h
	mex -largeArrayDims -g BetheApprox_mex.cpp MinSum.cpp BetheApprox.cpp graph.cpp maxflow.cpp -o BetheApprox_debug_mex.mexmaci64

clean:
	rm *.o *.mexmaci64
