emd_flow: main.cc emd_flow.o emd_flow.h emd_flow_network_factory.o emd_flow_network_factory.h emd_flow_network_sap.o
	g++ -Wall -Wextra -O2 -o emd_flow main.cc emd_flow.o emd_flow_network_factory.o emd_flow_network_sap.o -lboost_program_options -L lemon/lib -lemon

emd_flow.o: emd_flow.cc emd_flow.h emd_flow_network.h emd_flow_network_factory.h
	g++ -Wall -Wextra -fPIC -O2 -c -o emd_flow.o emd_flow.cc

emd_flow_network_factory.o: emd_flow_network_factory.cc emd_flow_network_factory.h emd_flow_network_lemon.h emd_flow_network_sap.h emd_flow_network.h
	g++ -Wall -Wextra -fPIC -O2 -c -o emd_flow_network_factory.o emd_flow_network_factory.cc -I lemon/include

emd_flow_network_sap.o: emd_flow_network_sap.cc emd_flow_network_sap.h emd_flow_network.h
	g++ -Wall -Wextra -fPIC -O2 -g -c -o emd_flow_network_sap.o emd_flow_network_sap.cc

mexfile: emd_flow.o emd_flow_network_factory.o emd_flow_network_sap.o emd_flow.h emd_flow_network_factory.h mex_wrapper.cc mex_helper.h
	mex -v CXXFLAGS="\$$CXXFLAGS -Wall -Wextra" -output emd_flow mex_wrapper.cc emd_flow.o emd_flow_network_sap.o emd_flow_network_factory.o

emd_flow_lambda_mexwrapper: emd_flow_network.o emd_flow_lambda_mexwrapper.cc
	mex -output emd_flow_lambda emd_flow_lambda_mexwrapper.cc emd_flow_network.o -Ilemon/include

all: emd_flow mexfile

clean:
	rm -f *.o
	rm -f emd_flow
	rm -f *.mexa64
