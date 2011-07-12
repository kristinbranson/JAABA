
#include <jsonrpc/jsonrpc.h>
#include <json/json.h>
#include <cstdio>
#include <cstdlib>
#include <csignal>
#include <assert.h>

#include <omp.h>

#include "../svm/svm_struct/svm_struct_api_behavior_sequence.h"
#include "../svm/svm_struct/online_structured_learning.h"
#include "svgPlotter.h"

#include "cv.h"
#include "highgui.h"
IplImage *VisualizeBouts(BehaviorBoutSequence *seq, BehaviorGroups *groups, int beh, const char *fname, char *html);


#define JSON_ERROR(errStr) \
  fprintf(stderr, "%s\n", errStr); \
  response["error"] = std::string(errStr); \
  return false;



class StructuredLearnerRpc
{
  StructuredSVMOnlineLearner *learner;

  omp_lock_t lock;

public:
  StructuredLearnerRpc();
  int main(int argc, const char **argv) ;

  bool AddNewExample(const Json::Value& root, Json::Value& response);
  bool ClassifyExample(const Json::Value& root, Json::Value& response);
  bool GetStatistics(const Json::Value& root, Json::Value& response);

  bool SaveCurrentModel(const Json::Value& root, Json::Value& response);
  bool Shutdown(const Json::Value& root, Json::Value& response);
  bool SetParameter(const Json::Value& root, Json::Value& response);

private:
  int RunServer(int port);

  void Lock() { omp_set_lock(&lock); }
  void Unlock() { omp_unset_lock(&lock); }

private:
};

StructuredLearnerRpc::StructuredLearnerRpc() {
  omp_init_lock(&lock);
}

bool StructuredLearnerRpc::Shutdown(const Json::Value& root, Json::Value& response) {
  if(learner) learner->Shutdown();
  
  return true;
}

bool StructuredLearnerRpc::SaveCurrentModel(const Json::Value& root, Json::Value& response) {
  char fname[1000]; strcpy(fname, root.get("filename", "").asString().c_str());
  if(learner) { 
    learner->SaveModel(strlen(fname) ? fname : NULL); 
    response["filename"] = fname;
  }
  return true;
}

bool StructuredLearnerRpc::SetParameter(const Json::Value& root, Json::Value& response) {
  double c = root.get("C", -1).asDouble();
  double lambda = root.get("lambda", -1).asDouble();
  if(lambda >= 0 && learner) { learner->SetLambda(lambda); response["lambda"] = lambda; }
  else if(c >= 0 && learner) { learner->SetC(c); response["C"] = c; }
  else { JSON_ERROR("SetParameter failed.  Supported parameters are 'C' or 'lambda'"); }
  return true;
}

bool StructuredLearnerRpc::GetStatistics(const Json::Value& root, Json::Value& response) {

  if(learner) {
    char test_train_by_example_plot[1000], error_decomp_by_example_plot[1000], test_train_by_iteration_plot[1000], error_decomp_by_iteration_plot[1000];
    strcpy(test_train_by_example_plot, root.get("test_train_by_example_plot_name", "").asString().c_str());
    strcpy(error_decomp_by_example_plot, root.get("error_decomp_by_example_plot_name", "").asString().c_str());
    strcpy(test_train_by_iteration_plot, root.get("test_train_by_iteration_plot_name", "").asString().c_str());
    strcpy(error_decomp_by_iteration_plot, root.get("error_decomp_by_iteration_plot_name", "").asString().c_str());

    int width = root.get("width", 0).asInt();
    int height = root.get("height", 0).asInt();
    const char *css = root.get("css", "").asString().c_str();
    int window = root.get("window", 100).asInt();
    char fname[1000], errStr[1000];

    if(strlen(test_train_by_example_plot)) {
      double *train_err_buff, *test_err_buff;
      SVGPlotter plotter;
      long n;
      if(strlen(css)) plotter.SetCSS(css);
      if(width && height) plotter.SetSize(width, height);
      learner->GetStatisticsByExample(window, &n, NULL, NULL, NULL, NULL, &train_err_buff, &test_err_buff);
      plotter.AddPlot(NULL, train_err_buff, n, 1, 1, "Train error", "class1");
      plotter.AddPlot(NULL, test_err_buff, n, 1, 1, "Test error", "class2");
      plotter.SetXLabel("Number of Training Examples (n)");
      plotter.SetYLabel("Error");
      plotter.SetTitle("Structured SVM Error");
      double x, y;
      if((x=root.get("xmin", -100000).asDouble()) != -100000) plotter.SetXMin(x);
      if((y=root.get("ymin", -100000).asDouble()) != -100000) plotter.SetYMin(y);
      if((x=root.get("xmax", -100000).asDouble()) != -100000) plotter.SetXMax(x);
      if((y=root.get("ymax", -100000).asDouble()) != -100000) plotter.SetYMax(y);
      sprintf(fname, "%s", test_train_by_example_plot);
	  if(!plotter.Save(fname)) { sprintf(errStr, "Error saving to %s\n", fname); JSON_ERROR(errStr); }
      StripFileExtension(fname); strcat(fname, ".m"); 
	  if(!plotter.Save(fname)) { sprintf(errStr, "Error saving to %s\n", fname); JSON_ERROR(errStr); }
      response["test_train_by_example_plot_name"] = fname;
      free(train_err_buff);  free(test_err_buff);
    }
    if(strlen(error_decomp_by_example_plot)) {
      double *gen_err_buff, *opt_err_buff, *model_err_buff, *reg_err_buff;
      SVGPlotter plotter;
      long n;
      if(strlen(css)) plotter.SetCSS(css);
      if(width && height) plotter.SetSize(width, height);
      learner->GetStatisticsByExample(window, &n, &gen_err_buff, &opt_err_buff, &model_err_buff, &reg_err_buff, NULL, NULL);
      plotter.AddPlot(NULL, gen_err_buff, n, 1, 1, "Test Error", "class1");
      plotter.AddPlot(NULL, opt_err_buff, n, 1, 1, "Optimization Error", "class2");
      plotter.AddPlot(NULL, model_err_buff, n, 1, 1, "Model Error", "class3");
      plotter.AddPlot(NULL, reg_err_buff, n, 1, 1, "Regularization Error", "class4");
      plotter.SetXLabel("Number of Training Examples (n)");
      plotter.SetYLabel("Error");
      plotter.SetTitle("Error Decomposition");
      double x, y;
      if((x=root.get("xmin", -100000).asDouble()) != -100000) plotter.SetXMin(x);
      if((y=root.get("ymin", -100000).asDouble()) != -100000) plotter.SetYMin(y);
      if((x=root.get("xmax", -100000).asDouble()) != -100000) plotter.SetXMax(x);
      if((y=root.get("ymax", -100000).asDouble()) != -100000) plotter.SetYMax(y);
      sprintf(fname, "%s", error_decomp_by_example_plot);
	  if(!plotter.Save(fname)) { sprintf(errStr, "Error saving to %s\n", fname); JSON_ERROR(errStr); }
      StripFileExtension(fname); strcat(fname, ".m"); 
	  if(!plotter.Save(fname)) { sprintf(errStr, "Error saving to %s\n", fname); JSON_ERROR(errStr); }
      response["error_decomp_by_example_plot_name"] = fname;
      free(gen_err_buff);  free(opt_err_buff); free(model_err_buff);  free(reg_err_buff); 
    }
    if(strlen(test_train_by_iteration_plot)) {
      double *train_err_buff, *test_err_buff, *time_buff;
      SVGPlotter plotter;
      long t, tm;
      if(strlen(css)) plotter.SetCSS(css);
      if(width && height) plotter.SetSize(width, height);
      learner->GetStatisticsByIteration(window, &t, &tm, NULL, NULL, NULL, NULL, &train_err_buff, &test_err_buff, &time_buff);
      plotter.AddPlot(time_buff, train_err_buff, t, 0, 0, "Train error", "class1");
      plotter.AddPlot(time_buff, test_err_buff, t, 0, 0, "Test error", "class2");
      plotter.SetXLabel("Training Time (seconds)");
      plotter.SetYLabel("Error");
      plotter.SetTitle("Structured SVM Error");
      sprintf(fname, "%s", test_train_by_iteration_plot);
      double x, y;
      if((x=root.get("xmin", -100000).asDouble()) != -100000) plotter.SetXMin(x);
      if((y=root.get("ymin", -100000).asDouble()) != -100000) plotter.SetYMin(y);
      if((x=root.get("xmax", -100000).asDouble()) != -100000) plotter.SetXMax(x);
      if((y=root.get("ymax", -100000).asDouble()) != -100000) plotter.SetYMax(y);
	  if(!plotter.Save(fname)) { sprintf(errStr, "Error saving to %s\n", fname); JSON_ERROR(errStr); }
      StripFileExtension(fname); strcat(fname, ".m"); 
	  if(!plotter.Save(fname)) { sprintf(errStr, "Error saving to %s\n", fname); JSON_ERROR(errStr); }
      response["test_train_by_iteration_plot_name"] = fname;
      free(train_err_buff);  free(test_err_buff); free(time_buff);
    }
    if(strlen(error_decomp_by_iteration_plot)) {
      double *gen_err_buff, *opt_err_buff, *model_err_buff, *reg_err_buff, *time_buff;
      SVGPlotter plotter;
      long t, tm;
      if(strlen(css)) plotter.SetCSS(css);
      if(width && height) plotter.SetSize(width, height);
      learner->GetStatisticsByIteration(window, &t, &tm, &gen_err_buff, &opt_err_buff, &model_err_buff, &reg_err_buff, NULL, NULL, &time_buff);
      plotter.AddPlot(time_buff, gen_err_buff, t, 0, 0, "Test Error", "class1");
      plotter.AddPlot(time_buff, opt_err_buff, t, 0, 0, "Optimization Error", "class2");
      plotter.AddPlot(time_buff, model_err_buff, t, 0, 0, "Model Error", "class3");
      plotter.AddPlot(time_buff, reg_err_buff, t, 0, 0, "Regularization Error", "class4");
      plotter.SetXLabel("Training Time (seconds)");
      plotter.SetYLabel("Error");
      plotter.SetTitle("Error Decomposition");
      double x, y;
      if((x=root.get("xmin", -100000).asDouble()) != -100000) plotter.SetXMin(x);
      if((y=root.get("ymin", -100000).asDouble()) != -100000) plotter.SetYMin(y);
      if((x=root.get("xmax", -100000).asDouble()) != -100000) plotter.SetXMax(x);
      if((y=root.get("ymax", -100000).asDouble()) != -100000) plotter.SetYMax(y);
      sprintf(fname, "%s", error_decomp_by_iteration_plot);
	  if(!plotter.Save(fname)) { sprintf(errStr, "Error saving to %s\n", fname); JSON_ERROR(errStr); }
      StripFileExtension(fname); strcat(fname, ".m"); 
	  if(!plotter.Save(fname)) { sprintf(errStr, "Error saving to %s\n", fname); JSON_ERROR(errStr); }
      response["error_decomp_by_iteration_plot_name"] = fname;
      free(gen_err_buff);  free(opt_err_buff); free(model_err_buff);  free(reg_err_buff); free(time_buff);
    }
  }

  return true;
}

bool StructuredLearnerRpc::ClassifyExample(const Json::Value& root, Json::Value& response) {
  if(learner) {
    char fname[1000]; strcpy(fname, root.get("input_name", "").asString().c_str());
    char oname[1000]; strcpy(oname, root.get("output_name", "").asString().c_str());
    char vname[1000]; strcpy(vname, root.get("visualization_name", "").asString().c_str());

	if(!strlen(fname)) {
      JSON_ERROR("No 'input_name' parameter passed to classify_example()\n");
	} else {
      SVMStructMethod *m = learner->GetStructMethod();
      STRUCTMODEL sm = *learner->GetStructModel();
      STRUCT_LEARN_PARM *sparm = learner->GetStructLearnParms();
      EXAMPLE ex = m->read_struct_example(fname, sparm);
	  BehaviorBoutSequence *bouts = (BehaviorBoutSequence*)ex.y.data;
	  BehaviorBoutFeatures *feat = (BehaviorBoutFeatures*)ex.x.data;
	  if(root.get("use_partial_label", true).asBool()) 
        feat->partial_label = bouts;
      sm.w = learner->GetCurrentWeights();
      double score;
      LABEL y = m->classify_struct_example(&ex.x, &sm, sparm, &score);
      BehaviorBoutSequence *preds = (BehaviorBoutSequence*)y.data;
	  if(strlen(oname)) 
        m->save_example(feat->data, preds, oname);
	  if(strlen(vname)) 
        VisualizeBouts(preds, preds->behaviors, 0, vname, NULL);
	  
      response["score"] = score;
	}
  }
  return true;
}

bool StructuredLearnerRpc::AddNewExample(const Json::Value& root, Json::Value& response) {
  if(learner) {
    char fname[1000]; strcpy(fname, root.get("filename", "").asString().c_str());
	if(!strlen(fname)) {
      JSON_ERROR("No 'fname' parameter passed to add_new_example()\n");
	} else {
	  learner->AddExample(fname);
	  response["filename"] = fname;

	}
  }
  return true;
}

/**
 * \var g_run
 * \brief Running state of the program.
 */
static volatile bool g_run = false;

/**
 * \brief Signal management.
 * \param code signal code
 */
static void signal_handler(int code)
{
  switch(code)
  {
    case SIGINT:
    case SIGTERM:
      g_run = false;
      break;
    default:
      break;
  }
}

int main(int argc, const char **argv) {
  StructuredLearnerRpc v;
  v.main(argc, argv);
}

int StructuredLearnerRpc::main(int argc, const char **argv) {
  STRUCT_LEARN_PARM struct_parm;
  STRUCTMODEL structmodel;

  memset(&structmodel, 0, sizeof(STRUCTMODEL));
  memset(&struct_parm, 0, sizeof(STRUCT_LEARN_PARM));
  struct_parm.method = SPO_DUAL_UPDATE_WITH_CACHE;

  //srand ( time(NULL) );

  /* to avoid compilation warnings */
  argc = argc;
  argv = argv;

  int port = 8086;
  
  
  int retval = -1;
  omp_set_nested(1);
  #pragma omp parallel num_threads(2) 
  {
    if(omp_get_thread_num() == 0)
      train_main(argc, argv, &struct_parm, &structmodel, NULL, &learner); 
    else 
      retval = RunServer(port);
  } 

  return retval;
}



int StructuredLearnerRpc::RunServer(int port) {
  Json::Rpc::TcpServer server(std::string("127.0.0.1"), port);

  if(!networking::init())
  {
    std::cerr << "Networking initialization failed" << std::endl;
    exit(EXIT_FAILURE);
  }

  if(signal(SIGTERM, signal_handler) == SIG_ERR)
  {
    std::cout << "Error signal SIGTERM will not be handled" << std::endl;
  }

  if(signal(SIGINT, signal_handler) == SIG_ERR)
  {
    std::cout << "Error signal SIGINT will not be handled" << std::endl;
  }

  server.AddMethod(new Json::Rpc::RpcMethod<StructuredLearnerRpc>(*this, &StructuredLearnerRpc::AddNewExample, std::string("add_example")));
  server.AddMethod(new Json::Rpc::RpcMethod<StructuredLearnerRpc>(*this, &StructuredLearnerRpc::ClassifyExample, std::string("classify_example")));
  server.AddMethod(new Json::Rpc::RpcMethod<StructuredLearnerRpc>(*this, &StructuredLearnerRpc::GetStatistics, std::string("plot_stats")));
  server.AddMethod(new Json::Rpc::RpcMethod<StructuredLearnerRpc>(*this, &StructuredLearnerRpc::SaveCurrentModel, std::string("save")));
  server.AddMethod(new Json::Rpc::RpcMethod<StructuredLearnerRpc>(*this, &StructuredLearnerRpc::Shutdown, std::string("shutdown")));
  server.AddMethod(new Json::Rpc::RpcMethod<StructuredLearnerRpc>(*this, &StructuredLearnerRpc::SetParameter, std::string("change_parameter")));
  


  fprintf(stderr, "Server listening on port %d\n", port);


  if(!server.Bind())
  {
    std::cout << "Bind failed" << std::endl;
    exit(EXIT_FAILURE);
  }

  if(!server.Listen())
  {
    std::cout << "Listen failed" << std::endl;
    exit(EXIT_FAILURE);
  }

  g_run = true;

  std::cout << "Start JSON-RPC TCP server" << std::endl;

  while(g_run)
  {
    server.WaitMessage(1000);
  }

  std::cout << "Stop JSON-RPC TCP server" << std::endl;
  server.Close();
  networking::cleanup();

  return EXIT_SUCCESS;
  
}


