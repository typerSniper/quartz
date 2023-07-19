#include "quartz/parser/qasm_parser.h"
#include "quartz/tasograph/substitution.h"
#include "quartz/tasograph/tasograph.h"

using namespace quartz;

void parse_args(char **argv, int argc, bool &simulated_annealing,
                bool &early_stop, bool &disable_search, int& timeout,
                std::string &input_filename, std::string &output_filename,
                std::string &eqset_filename) {
  assert(argv[1] != nullptr);
  input_filename = std::string(argv[1]);
  early_stop = true;
  for (int i = 2; i < argc; i++) {
    if (!std::strcmp(argv[i], "--output")) {
      output_filename = std::string(argv[++i]);
      continue;
    }
    if (!std::strcmp(argv[i], "--eqset")) {
      eqset_filename = std::string(argv[++i]);
      continue;
    }
    if (!std::strcmp(argv[i], "--disable_search")) {
      disable_search = true;
      continue;
    }
    if (!std::strcmp(argv[i], "--timeout")) {
      timeout = atoi(argv[++i]);
      continue;
    }
    if (!std::strcmp(argv[i], "--outfile")) {
      output_filename = std::string(argv[++i]);
      continue;
    }
  }
}

auto clifft_set = {GateType::t, GateType::tdg, GateType::h, GateType::s,
               GateType::sdg, GateType::cx,  GateType::input_qubit};
auto tcount = [](Graph *graph) {
    return (float) (graph->specific_gate_count(GateType::t) + graph->specific_gate_count(GateType::tdg));};

auto cost = [](Graph *graph) {
    return (float) (graph->total_cost() + 0.0 * graph->specific_gate_count(GateType::t) + 0.0 * graph->specific_gate_count(GateType::tdg));};


int main(int argc, char **argv) {
  std::string input_fn, output_fn;
  std::string eqset_fn = "../Nam_6_3_complete_ECC_set.json";
  bool simulated_annealing = false;
  bool early_stop = false;
  bool disable_search = false;
  int timeout = 60;
  parse_args(argv, argc, simulated_annealing, early_stop, disable_search, timeout,
             input_fn, output_fn, eqset_fn);
  auto fn = input_fn.substr(input_fn.rfind('/') + 1);

  // Construct contexts
  Context src_ctx(clifft_set);
  // Load qasm file
  QASMParser qasm_parser(&src_ctx);
  CircuitSeq *dag = nullptr;
  if (!qasm_parser.load_qasm(input_fn, dag)) {
    std::cout << "Parser failed" << std::endl;
  }
  Graph graph(&src_ctx, dag);
  auto graph_before_search = &graph;

  auto start = std::chrono::steady_clock::now();

  std::shared_ptr<Graph> graph_after_search;
  if (disable_search) {
    std::cout << "calling gopt" <<std::endl;
    graph_after_search = graph_before_search->greedy_optimize(&src_ctx, eqset_fn, true, cost, timeout);
  }
  else {
    std::cout << "general search" << std::endl;
    auto t0 = std::chrono::steady_clock::now();
    auto greedy_graph = graph_before_search->greedy_optimize(&src_ctx, eqset_fn, true, cost, timeout);
    auto t1 = std::chrono::steady_clock::now();
    auto rem = timeout - std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count()/1000.0;
    std::cout << "running search with remaining time " << rem << std::endl;
    if (rem > 0) {
      graph_after_search = greedy_graph->optimize(&src_ctx, eqset_fn, fn, true, cost, -1, rem);
    }
    else {
      graph_after_search = greedy_graph;
      std::cout << "(0.0, " << cost (graph_after_search.get()) << ");\n";
    }
  }
  auto end = std::chrono::steady_clock::now();
  std::cout << "Optimization results of Quartz for " << fn
            << " on Clifford T gate set."
            << " Cost after optimization: "
            << cost(graph_after_search.get()) << ", "
            << "T count: " << tcount(graph_after_search.get()) << ", "
            << "before T count: " << tcount(graph_before_search) << ", "
            << (double)std::chrono::duration_cast<std::chrono::milliseconds>(
                   end - start)
                       .count() /
                   1000.0
            << " seconds." << std::endl;
  graph_after_search->to_qasm(output_fn, false, false);
}
