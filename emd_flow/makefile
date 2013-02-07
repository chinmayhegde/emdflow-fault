emd_flow: main.cc emd_flow.o emd_flow.h emd_flow_network.o
	g++ -Wall -Wextra -O2 -o emd_flow main.cc emd_flow.o emd_flow_network.o -lboost_program_options

emd_flow.o: emd_flow.cc emd_flow.h emd_flow_network.h
	g++ -Wall -Wextra -fPIC -O2 -c -o emd_flow.o emd_flow.cc -I lemon/include

emd_flow_network.o: emd_flow_network.cc emd_flow_network.h
	g++ -Wall -Wextra -fPIC -O2 -c -o emd_flow_network.o emd_flow_network.cc -lemon -I lemon/include -L lemon/lib

mexfile: emd_flow.o emd_flow_network.o emd_flow.h mex_wrapper.cc
	mex -output emd_flow mex_wrapper.cc emd_flow.o emd_flow_network.o

emd_flow_lambda_mexwrapper: emd_flow_network.o emd_flow_lambda_mexwrapper.cc
	mex -output emd_flow_lambda emd_flow_lambda_mexwrapper.cc emd_flow_network.o -Ilemon/include

all: emd_flow mexfile

clean:
	rm -f *.o
	rm -f emd_flow
	rm -f *.mexa64
