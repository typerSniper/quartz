#include "quartz/parser/qasm_parser.h"
#include "quartz/tasograph/substitution.h"
#include "quartz/tasograph/tasograph.h"
#include "quartz_api.h"
using namespace quartz;

// Construct contexts
auto nam_set = {GateType::h, GateType::x, GateType::rz, GateType::add,
                  GateType::cx, GateType::input_qubit, GateType::input_param};
auto clifft_set = {GateType::t, GateType::tdg, GateType::h, GateType::s,
               GateType::sdg, GateType::cx,  GateType::input_qubit};
auto gate_set = nam_set;
// Context gctxt();


// Context ctctxt({GateType::t, GateType::tdg, GateType::h, GateType::s, GateType::sdg, GateType::cx});

// Context* new_ctxt (const Context c) {
//   auto vec = c.get_supported_gates();
//   return new Context (vec);
// }

auto total_count = [](Graph *graph) { return graph->total_cost(); };
auto tcount = [](Graph *graph) {
    return (float) (graph->specific_gate_count(GateType::t) + graph->specific_gate_count(GateType::tdg));};

auto gcost_function = total_count;

int write_qasm_to_buffer (std::string cqasm, char* buffer, int buff_size) {
  int blen = static_cast<int>(strlen(cqasm.c_str()));
  if (blen > buff_size) {
    return -1 * blen;
  } else {
    strcpy(buffer, cqasm.c_str());
    return blen;
  }
}

extern "C" long unsigned int load_eqset_ (const char* eqset_fn_, unsigned char** store) {
  std::string eqset_fn(eqset_fn_);
  EquivalenceSet* eqs = new EquivalenceSet();
  ParamInfo* param_info = new ParamInfo();
  auto gctxt = Context(gate_set, param_info);
  if (!eqs->load_json(&gctxt, eqset_fn, false)) {
    std::cout << "Failed to load equivalence file \"" << eqset_fn
              << "\"." << std::endl;
    assert(false);
  }
  std::vector<std::vector<CircuitSeq *>> eccs = eqs->get_all_equivalence_sets();
  std::vector<CircuitSeq *> *arr = new std::vector<CircuitSeq*>[eccs.size()];
  std::copy(eccs.begin(), eccs.end(), arr);
  *store = reinterpret_cast<unsigned char*> (arr);
  return eccs.size();
}

std::vector<GraphXfer *> filter_greedy (Context * ctxt, std::vector<std::vector<CircuitSeq *>>& eccs) {
  std::vector<GraphXfer *> xfers;
  for (const auto &ecc : eccs) {
    const int ecc_size = (int)ecc.size();
    std::vector<Graph> graphs;
    std::vector<int> graph_cost;
    graphs.reserve(ecc_size);
    graph_cost.reserve(ecc_size);
    for (auto &circuit : ecc) {
      graphs.emplace_back(ctxt, circuit);
      graph_cost.emplace_back(gcost_function(&(graphs.back())));
    }
    int representative_id =
        (int)(std::min_element(graph_cost.begin(), graph_cost.end()) -
              graph_cost.begin());
    for (int i = 0; i < ecc_size; i++) {
      if (graph_cost[i] != graph_cost[representative_id]) {
        auto xfer = GraphXfer::create_GraphXfer(ctxt, ecc[i],
                                                ecc[representative_id], true);
        if (xfer != nullptr) {
          xfers.push_back(xfer);
        }
      }
    }
  }
  return xfers;
}

extern "C" long unsigned int load_greedy_xfers_ (const char* eqset_fn_, unsigned char** store) {
  // std::string eqset_fn(eqset_fn_);
  std::string eqset_fn(reinterpret_cast<const char*>(eqset_fn_));
  std::istringstream iss(eqset_fn);

  EquivalenceSet* eqs = new EquivalenceSet();
  ParamInfo* param_info = new ParamInfo();
  Context *ctxt = new Context(gate_set, param_info);
  if (!eqs->load_json(ctxt, iss, false)) {
    std::cout << "Failed to load equivalence file \"" << "whe"
              << "\"." << std::endl;
    assert(false);
  }
  std::vector<std::vector<CircuitSeq *>> eccs = eqs->get_all_equivalence_sets();
  std::vector<GraphXfer *> xfers = filter_greedy(ctxt, eccs);

  // std::cout << "greedy_optimize(): Number of xfers that reduce cost: "
              // << xfers.size() << std::endl;

  std::vector<GraphXfer *> *vptr = new std::vector<GraphXfer*>(xfers);
  *store = reinterpret_cast<unsigned char*> (vptr);
  return xfers.size();
}

extern "C" void load_xfers_ (const char* eqset_fn_,
  unsigned char** gstore, long unsigned int* glen,
  unsigned char** astore, long unsigned int* alen) {
  // std::string eqset_fn(eqset_fn_);
  std::string eqset_fn(reinterpret_cast<const char*>(eqset_fn_));
  std::istringstream iss(eqset_fn);

  EquivalenceSet* eqs = new EquivalenceSet();
  ParamInfo* param_info = new ParamInfo();
  Context *ctxt = new Context (gate_set, param_info);
  if (!eqs->load_json(ctxt, iss, false)) {
    std::cout << "Failed to load equivalence file \"" << "whe"
              << "\"." << std::endl;
    assert(false);
  }
  std::vector<std::vector<CircuitSeq *>> eccs = eqs->get_all_equivalence_sets();
  std::vector<GraphXfer *> xfers;

  for (const auto &ecc : eccs) {
    CircuitSeq *representative = ecc.front();
    for (auto &circuit : ecc) {
      if (circuit != representative) {
        auto xfer =
            GraphXfer::create_GraphXfer(ctxt, circuit, representative, true);
        if (xfer != nullptr) {
          xfers.push_back(xfer);
        }
        xfer = GraphXfer::create_GraphXfer(ctxt, representative, circuit, true);
        if (xfer != nullptr) {
          xfers.push_back(xfer);
        }
      }
    }
  }
  std::cout << "Number of xfers: " << xfers.size() << std::endl;
  std::vector<GraphXfer *> *aptr = new std::vector<GraphXfer*>(xfers);
  *astore = reinterpret_cast<unsigned char*> (aptr);
  *alen =xfers.size();

  auto greedy_xfers = filter_greedy (ctxt, eccs);
  std::vector<GraphXfer *> *gptr = new std::vector<GraphXfer*>(greedy_xfers);
  *gstore = reinterpret_cast<unsigned char*> (gptr);
  *glen = greedy_xfers.size();

  return;
}

extern "C" int preprocess_ (const char* cqasm_, char* buffer, int buff_size) {

  std::string cqasm(cqasm_);
  ParamInfo* param_info = new ParamInfo();
  Context src_ctx({GateType::h, GateType::ccz, GateType::x, GateType::cx, GateType::rz,
                   GateType::input_qubit, GateType::input_param}, param_info);



  QASMParser qasm_parser(&src_ctx);
  CircuitSeq *dag = nullptr;
  if (!qasm_parser.load_qasm_str(cqasm, dag)) {
    std::cout << "Parser failed" << std::endl;
    return -1;
  }
  auto graph = Graph::from_qasm_str (&src_ctx, cqasm);
  ParamInfo* param_info2 = new ParamInfo();
  Context dst_ctx({GateType::h, GateType::x, GateType::rz, GateType::add,
                   GateType::cx, GateType::input_qubit, GateType::input_param}, param_info2);

  auto union_ctx = union_contexts(&src_ctx, &dst_ctx);
  auto xfer_pair = GraphXfer::ccz_cx_rz_xfer(&src_ctx, &dst_ctx, &union_ctx);
  auto new_graph = graph->toffoli_flip_greedy(GateType::rz, xfer_pair.first, xfer_pair.second);
  new_graph->constant_and_rotation_elimination();

  std::string new_qasm = new_graph->to_qasm(false, false);
  return write_qasm_to_buffer (new_qasm, buffer, buff_size);

  // decompose ccz as cx and rz
  // Context rem_ctx({GateType::u1, GateType::rx, GateType::h, GateType::x, GateType::rz, GateType::add,
  //                  GateType::cx, GateType::input_qubit, GateType::input_param});
  // auto imt_ctx = union_contexts(&src_ctx, &rem_ctx);
  // // std::cout << "flipping done\n"<<  std::endl;

  // Context dst_ctx(gate_set);
  // auto uctx = union_contexts(&rem_ctx, &dst_ctx);
  // // auto y_rule = "y = rz(0.5pi) q0; rz(0.5pi) q0; h q0; rz(0.5pi) q0; rz(0.5pi) q0; h q0;";
  // RuleParser rules({"rx q0 p0 = h q0; rz q0 p0; h q0;", "u1 q0 p0 = rz q0 p0;"}); // TODO: check this.
  // // std::cout << "contexy shifting " << new_graph->to_qasm (false, false) <<  std::endl;
  // auto fin_graph = new_graph->context_shift(&rem_ctx, &dst_ctx, &uctx, &rules, false);

  // std::cout << "ruling done\n"<<  std::endl;

}


extern "C" int opt_circuit_ (const char* cqasm_, int timeout, char* buffer, int buff_size, unsigned char* xfers_) {

  std::string cqasm(cqasm_);

  std::vector<GraphXfer*>* xfers_ptr = reinterpret_cast<std::vector<GraphXfer*>*> (xfers_);
  std::vector<GraphXfer*> xfers = *xfers_ptr;

  // std::string eqset_fn = "Nam_4_3_complete_ECC_set.json";
  Context* ctxt;
  if (xfers.size () == 0) {
    // ctxt = new_ctxt(gctxt) ;
    ParamInfo* param_info = new ParamInfo();
    ctxt = new Context (gate_set, param_info);
  } else {
    std::cout << "reusing context\n";
    ctxt = xfers[0]->dst_ctx_;
  }
  std::cout << "pinfo " << ctxt->param_info_to_json() << std::endl;
  auto graph = Graph::from_qasm_str (ctxt, cqasm);

  auto start = std::chrono::steady_clock::now();
  // Assume that the context is same?
  // std::cout << "calling greedy_opt" << std::endl;
  // auto graph_after_search = graph.greedy_optimize(ctxt, eqset_fn, /*print_message=*/ false);
  std::cout << "timeout received = " << timeout << std::endl;
  std::shared_ptr<Graph> graph_after_search;
  std::shared_ptr<Graph> graph_after_search2;
  timeout = 0;
  if (timeout == 0) {
    std::cout << "calling here\n";
    graph_after_search = graph->greedy_optimize(ctxt, "lib/quartz/Nam_6_3_complete_ECC_set.json", /*print_message=*/ false, gcost_function);
    std::cout << "pinfo " << ctxt->param_info_to_json() << std::endl;

    graph_after_search = graph->greedy_optimize_with_xfers(ctxt, xfers, /*print_message=*/ false, gcost_function);
  }
  else if (timeout < 0) {
    graph_after_search = graph->greedy_optimize_with_xfers(ctxt, xfers, /*print_message=*/ false, gcost_function, -1 * timeout);
    std::cout << "done" << std::endl;
  }
  else {
    graph_after_search = graph->optimize(xfers, 1.05 * gcost_function(graph.get()), "", "", false, nullptr, timeout);
  }

  auto end = std::chrono::steady_clock::now();

  std::cout << " Cost function optimized from: "
            << gcost_function (graph.get()) << " to "
            << gcost_function(graph_after_search.get()) << ", "
            // << "third = " << gcost_function(graph_after_search2.get()) << " "
            << (double)std::chrono::duration_cast<std::chrono::milliseconds>(
                   end - start)
                       .count() /
                   1000.0
            << " seconds." << std::endl;

  if (gcost_function(graph.get()) <= gcost_function(graph_after_search.get())) {
    return -1;
  }
  std::string cqasm2 = graph_after_search->to_qasm(false, false);

  if (xfers.size() == 0) {
    std::cout << "deleting wrongs" << std::endl;
    delete ctxt;
  }

  *xfers_ptr = xfers;
  exit(1);
  return write_qasm_to_buffer(cqasm2, buffer, buff_size);
  // std::cout << "circuit after opt = ";
  // std::cout << cqasm2.c_str() << std::endl;
  // int blen = static_cast<int>(strlen(cqasm2.c_str()));
  // if (blen > buff_size) {
  //   return -1 * blen;
  // } else {
  //   strcpy(buffer, cqasm2.c_str());
  //   return blen;
  // }
}
