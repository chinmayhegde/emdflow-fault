#include "emd_flow_network_factory.h"
#include "emd_flow_network.h"
#include "emd_flow_network_lemon.h"
#include "emd_flow_network_sap.h"

#include <memory>

#include <lemon/network_simplex.h>
#include <lemon/cost_scaling.h>
#include <lemon/capacity_scaling.h>

using namespace std;
using namespace lemon;

auto_ptr<EMDFlowNetwork> EMDFlowNetworkFactory::create_EMD_flow_network(
        const vector<vector<double> >& amplitudes, EMDFlowNetworkType type) {
  if (type == kLemonCostScaling) {
    return auto_ptr<EMDFlowNetwork>(
        new EMDFlowNetworkLemon<CostScaling<ListDigraph, int, double> >(
            amplitudes));
  } else if (type == kLemonNetworkSimplex) {
    return auto_ptr<EMDFlowNetwork>(
        new EMDFlowNetworkLemon<NetworkSimplex<ListDigraph, int, double> >(
            amplitudes));
  } else if (type == kLemonCapacityScaling) {
    return auto_ptr<EMDFlowNetwork>(
        new EMDFlowNetworkLemon<CapacityScaling<ListDigraph, int, double> >(
            amplitudes));
  } else if (type == kShortestAugmentingPath) {
    return auto_ptr<EMDFlowNetwork>(new EMDFlowNetworkSAP(amplitudes));
  } else {
    return auto_ptr<EMDFlowNetwork>();
  }
}

EMDFlowNetworkFactory::EMDFlowNetworkType EMDFlowNetworkFactory::parse_type(
    const std::string& name) {
  if (name == "lemon-costscaling") {
    return kLemonCostScaling;
  } else if (name == "lemon-networksimplex") {
    return kLemonNetworkSimplex;
  } else if (name == "lemon-capacityscaling") {
    return kLemonCapacityScaling;
  } else if (name == "sap" || name == "shortest-augmenting-path") {
    return kShortestAugmentingPath;
  } else {
    return kUnknownType;
  }
}


