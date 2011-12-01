#include "online_interactive_server.h"

#define SESSION_DIR "sessions"



bool StructuredLearnerRpc::FindOrCreateSession(const Json::Value& root, Json::Value& response) {
  int sess_ind = -1;
  char sess_id[1000];
  char errStr[1000];
  if(strlen(root.get("session_id", "").asString().c_str()) < 1000)
    strcpy(sess_id, root.get("session_id", "").asString().c_str());
  else
    strcpy(sess_id, "");
  
  if(!strlen(sess_id)) {
    if(!NewSession(root, response) || !InitializeSession(root, response)) {
      JSON_ERROR("Failed to create new session in classify_example\n", -1);
    }
    sess_ind = FindSession(response.get("session_id", "").asString().c_str());
    assert(sess_ind >= 0);
  } else {
    sess_ind=FindSession(sess_id, true);
    if(sess_ind < 0) {
      sprintf(errStr, "Invalid session_id %s in classify_example()\n", sess_id);  JSON_ERROR(errStr, sess_ind); 
    }
  }
  response["session_id"] = sessions[sess_ind].id;
  return sess_ind;
}

bool StructuredLearnerRpc::ClassifyExample(const Json::Value& root, Json::Value& response) {
  if(learner) {
    int sess_ind = FindOrCreateSession(root, response);
    if(sess_ind < 0) 
      return false;
    
    if(root.isMember("partial_label")) {
      if(!sessions[sess_ind].partial_label)
        sessions[sess_ind].partial_label = learner->NewStructuredLabel(sessions[sess_ind].example->x);
      if(!sessions[sess_ind].partial_label->load(root["partial_label"], learner)) {
	JSON_ERROR("Invalid 'partial_label' parameter", sess_ind); 
      }
    }

    SparseVector *w = learner->GetCurrentWeights();
    double score = learner->Inference(sessions[sess_ind].example->x, sessions[sess_ind].example->y, w, 
                               sessions[sess_ind].partial_label);
    
    if(!isnan(score)) response["score"] = score;
    response["y"] = sessions[sess_ind].example->y->save(learner);
    delete w;

    UnlockSession(sess_ind);

    return true;
  } else
    return false;
}

bool StructuredLearnerRpc::AddNewExample(const Json::Value& root, Json::Value& response) {
  if(learner) {
    char sess_id[1000];
    int sess_ind;
    if(strlen(root.get("session_id", "").asString().c_str()) < 1000)
      strcpy(sess_id, root.get("session_id", "").asString().c_str());
    else
      strcpy(sess_id, "");

    StructuredExample *ex = NULL;
    sess_ind = FindSession(sess_id, true);
    if(sess_ind >= 0) {
      ex = sessions[sess_ind].example;
      response["session_id"] = sessions[sess_ind].id;
    } else {
      ex = new StructuredExample;
      ex->x = learner->NewStructuredData();
      ex->y = learner->NewStructuredLabel(ex->x);
    }

    if(root.isMember("x")) {
      if(!ex->x->load(root["x"], learner)) { JSON_ERROR("Invalid 'x' parameter", sess_ind); }
    }
    if(root.isMember("y")) {
      if(!ex->y->load(root["y"], learner)) { JSON_ERROR("Invalid 'x' parameter", sess_ind); }
    } else if(sess_ind < 0) {
      delete ex;
      JSON_ERROR("No 'y' parameter specified", sess_ind); 
    }

    int ind = learner->AddExample(ex->x, ex->y);
    if(ind < 0) { JSON_ERROR("Error adding example\n", sess_ind) }
    response["index"] = ind;

    learner->SaveTrainingSet(ind);

    if(sess_ind >= 0)
      UnlockSession(sess_ind);
    else {
      delete ex;
    }

    return true;
  } else
    return false;
}





StructuredLearnerRpc::StructuredLearnerRpc(StructuredSVM *l) {
  num_sessions = 0;
  timestamp = 0;
  omp_init_lock(&lock);
  learner = l;
  server = NULL;
  train = true;
  runServer = false;
}

StructuredLearnerRpc::~StructuredLearnerRpc() {
  while(num_sessions) {
    EvictSession(0);
  }
  if(learner) {
    learner->Lock();
    delete learner;
  }
  if(server) delete server;
}

int StructuredLearnerRpc::FindSession(const char *sess_id, bool lock) {
  Lock();

  for(int i = 0; i < num_sessions; i++) {
    if(!strcmp(sessions[i].id, sess_id)) {
      sessions[i].timestamp = timestamp++;
      Unlock();
      if(lock) LockSession(i);
      return i;
    }
  }
  Unlock();
  return -1;
}

int StructuredLearnerRpc::FindSession(const Json::Value& root, Json::Value& response) {
  // Find the current session handle
  char sess_id[1000];
  if(strlen(root.get("session_id", "").asString().c_str()) < 1000)
    strcpy(sess_id, root.get("session_id", "").asString().c_str());
  else 
    return -1;
  
  int sess_ind=FindSession(sess_id, true);
  if(sess_ind < 0) 
    return -1;
  response["session_id"] = sessions[sess_ind].id;
  return sess_ind;
}

void StructuredLearnerRpc::EvictSession(int ind) {
  Lock();
  if(ind == -1) {
    unsigned long b = timestamp+1;
    for(int i = 0; i < num_sessions; i++) {
      if(sessions[i].timestamp < b) {
	b = sessions[i].timestamp;
	ind = i;
      }
    }
    assert(ind >= 0);
  }
  Unlock();

  LockSession(ind);
  delete sessions[ind].example;
  if(sessions[ind].partial_label)
    delete sessions[ind].partial_label;
  free(sessions[ind].id);
  omp_unset_lock(&sessions[ind].lock);

  Lock();
  for(int i = ind; i < num_sessions-1; i++)
    sessions[i] = sessions[i+1];
  num_sessions--;
  Unlock();
}


bool StructuredLearnerRpc::NewSession(const Json::Value& root, Json::Value& response) {
  // Randomly create a new session id
  Session s;
  memset(&s, 0, sizeof(Session));
  s.timestamp = timestamp++;
  char session_id[1000];
  do {
    sprintf(session_id, "%d", rand());
  } while(FindSession(session_id) >= 0);

  char sessDir[1000];
  sprintf(sessDir, "%s/%s", session_dir, session_id);

  if(root.get("mkdir", false).asBool()) 
     CreateDirectoryIfNecessary(sessDir, 777);
  
  s.id = StringCopy(session_id);
  s.example = new StructuredExample;
  s.example->x = learner->NewStructuredData();
  s.example->y = learner->NewStructuredLabel(s.example->x);

  if(num_sessions >= MAX_SESSIONS)
    EvictSession();
  omp_init_lock(&s.lock);
  sessions[num_sessions++] = s;

  response["session_id"] = session_id;
  response["session_dir"] = sessDir;

  return true;
}

bool StructuredLearnerRpc::InitializeSession(const Json::Value& root, Json::Value& response) {
  char errStr[1000];

  int sess_ind = -1;
  char sess_id[1000];
  strcpy(sess_id, "");
  if(strlen(root.get("session_id", "").asString().c_str()) < 1000)
    strcpy(sess_id, root.get("session_id", "").asString().c_str());
  if(!strlen(sess_id) && strlen(root.get("response", "").asString().c_str()) < 1000)
    strcpy(sess_id, response.get("session_id", "").asString().c_str());
  
  if((sess_ind=FindSession(sess_id, true)) < 0) { 
    sprintf(errStr, "Invalid session_id %s in InitializeSession()\n", sess_id);  JSON_ERROR(errStr, sess_ind); 
  } else {
    response["session_id"] = sessions[sess_ind].id;
    if(!sessions[sess_ind].example->x->load(root["x"], learner)) { JSON_ERROR("Invalid 'x' parameter", sess_ind); }
  }
  
  UnlockSession(sess_ind);

  return true;
}

bool StructuredLearnerRpc::Shutdown(const Json::Value& root, Json::Value& response) {
  if(learner) learner->Shutdown();
  response["status"] = "Shutting down";
  
  return true;
}

bool StructuredLearnerRpc::SaveCurrentModel(const Json::Value& root, Json::Value& response) {
  if(!IsSafeFileName(root.get("filename", "").asString().c_str())) { JSON_ERROR("Invalid 'filename' parameter to SaveCurrentModel()", -1); return false; }  
  char fname[1000]; strcpy(fname, root.get("filename", "").asString().c_str());
  if(learner) {
    bool saveFull = root.get("saveFull", true).asBool();
    response["filename"] = learner->Save(strlen(fname) ? fname : NULL, saveFull);
    return true;
  } else
    return false;
}

bool StructuredLearnerRpc::GetStatistics(const Json::Value& root, Json::Value& response) {
  /*if(learner) {
    char test_train_by_example_plot[1000], error_decomp_by_example_plot[1000], test_train_by_iteration_plot[1000], error_decomp_by_iteration_plot[1000];
    strcpy(test_train_by_example_plot, root.get("test_train_by_example_plot_name", "").asString().c_str());
    strcpy(error_decomp_by_example_plot, root.get("error_decomp_by_example_plot_name", "").asString().c_str());
    strcpy(test_train_by_iteration_plot, root.get("test_train_by_iteration_plot_name", "").asString().c_str());
    strcpy(error_decomp_by_iteration_plot, root.get("error_decomp_by_iteration_plot_name", "").asString().c_str());

    int width = root.get("width", 0).asInt();
    int height = root.get("height", 0).asInt();
    const char *css = root.get("css", "").asString().c_str();
    int window = root.get("window", 200).asInt();
    char fname[1000];

    if(strlen(test_train_by_example_plot) && IsSafeFileName(test_train_by_example_plot)) {
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
      sprintf(fname, "%s/%s", session_dir, test_train_by_example_plot);
      plotter.Save(fname);
      StripFileExtension(fname); strcat(fname, ".m"); plotter.Save(fname);
      response["test_train_by_example_plot_name"] = fname;
      free(train_err_buff);  free(test_err_buff);
    }
    if(strlen(error_decomp_by_iteration_plot) && IsSafeFileName(error_decomp_by_iteration_plot)) {
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
      sprintf(fname, "%s/%s", session_dir, error_decomp_by_iteration_plot);
      plotter.Save(fname);
      StripFileExtension(fname); strcat(fname, ".m"); plotter.Save(fname);
      response["error_decomp_by_iteration_plot_name"] = fname;
      free(gen_err_buff);  free(opt_err_buff); free(model_err_buff);  free(reg_err_buff); free(time_buff);
    }
  }
  */

  return true;
}


bool StructuredLearnerRpc::EvaluateTestset(const Json::Value& root, Json::Value& response) {
  
  if(!IsSafeFileName(root.get("testset", "").asString().c_str()))
     JSON_ERROR("Invalid 'testset' parameter passed to evaluate_testset()\n", -1);
  if(!IsSafeFileName(root.get("predictions", "").asString().c_str()))
     JSON_ERROR("Invalid 'predictions' parameter passed to evaluate_testset()\n", -1);

  char testset[1000];  strcpy(testset, root.get("testset", "").asString().c_str());
  char predictions[1000];  strcpy(predictions, root.get("predictions", "").asString().c_str());

  if(!strlen(testset)) {
    JSON_ERROR("No 'testset' parameter passed to evaluate_testset()\n", -1);
  } else if(learner) {
    response["ave_loss"] = learner->Test(testset, predictions);
  }
  return true;
}

bool StructuredLearnerRpc::SetParameter(const Json::Value& root, Json::Value& response) {
  double c = root.get("C", -1).asDouble();
  double lambda = root.get("lambda", -1).asDouble();
  double featureScale = root.get("feature_scale", 0).asDouble();
  int num_iter = root.get("num_iter", 5).asInt();
  bool set = false; 

  if(lambda >= 0 && learner) { learner->SetLambda(lambda, num_iter); response["lambda"] = lambda; set = true;  }
  else if(c >= 0 && learner) { learner->SetC(c, num_iter); response["C"] = c; set = true;  }

  if(featureScale > 0) { learner->SetFeatureScale(c, num_iter); response["feature_scale"] = featureScale; set = true; }

  if(!set) { JSON_ERROR("SetParameter failed.  Supported parameters are 'C' or 'lambda'", -1); }
  return true;
}





void StructuredLearnerRpc::parse_command_line_arguments(int argc, const char **argv) {
  sprintf(session_dir, SESSION_DIR);
  
  strcpy(testfile, "");
  strcpy(trainfile, "");
  strcpy(infile, "");
  strcpy(outfile, "");
  strcpy(paramfile, "");
  strcpy(predictionsfile, "");
  port = 8086;

  int i = 1;
  while(i < argc) {
    if(!strcmp(argv[i], "-P")) { 
      // Manually specify the port number
      port = atoi(argv[i+1]);
      runServer = true;
      i += 2;
    } else if(!strcmp(argv[i], "-s")) { 
      // Manually specify the directory where temporary files are stored
      strcpy(session_dir, argv[i+1]);
      i += 2;
    } else if(!strcmp(argv[i], "-p")) {
      // Read parameters for training
      assert(i+2 <= argc);
      strcpy(paramfile, argv[i+1]);
      i += 2;
    } else if(!strcmp(argv[i], "-d")) {
      // Optionally, start with an initial training set
      assert(i+2 <= argc);
      strcpy(trainfile, argv[i+1]);
      i += 2;
    } else if(!strcmp(argv[i], "-t")) {
      // Evaluate performance on a testset
      assert(i+2 <= argc);
      strcpy(testfile, argv[i+1]);
      if(i+3 <= argc && argv[i+2][0] != '-') { strcpy(predictionsfile, argv[i+2]); i++; }
      i += 2;
    } else if(!strcmp(argv[i], "-i")) {
      // Continue running from an in progress online learner that was previously saved to disk
      assert(i+2 <= argc);
      strcpy(infile, argv[i+1]);
      i += 2;
    } else if(!strcmp(argv[i], "-C")) {
      // Regularization parameter
      assert(i+2 <= argc);
      learner->SetC(atof(argv[i+1]));
      i += 2;
    } else if(!strcmp(argv[i], "-e")) {
      // Regularization parameter
      assert(i+2 <= argc);
      learner->SetEpsilon(atof(argv[i+1]));
      i += 2;
    } else if(!strcmp(argv[i], "-o")) {
      // Optionally, specify the filename to store the learned model
      assert(i+2 <= argc);
      strcpy(outfile, argv[i+1]);
      i += 2;
    } else if(argv[i][0] == '-') {
      i += 2;
    } else
      i++;
    /*else {
      fprintf(stderr, "Too many arguments\n"); 
      return -1; 
    } */
  }
  if(runServer) server = new JsonRpcServer(port);
}

int StructuredLearnerRpc::main(int argc, const  char **argv) {
  //srand ( time(NULL) );

  /* to avoid compilation warnings */
  argc = argc;
  argv = argv; 

  parse_command_line_arguments(argc, argv);

  if(strlen(infile)) {
    bool b = learner->Load(infile, true);
    assert(b);
  } else if(strlen(paramfile)) {
    bool b = learner->Load(paramfile, false);
    assert(b);
  } 

  if(strlen(trainfile)) {
    learner->LoadTrainset(trainfile);
    learner->GetTrainset()->Randomize();
  }

  int retval = -1;
  if(!runServer) {
    if(strlen(trainfile))   
      learner->Train(strlen(outfile) ? outfile : NULL, false);
    if(strlen(testfile)) 
      learner->Test(testfile, strlen(predictionsfile) ? predictionsfile : NULL);
  } else if(train) {
    omp_set_nested(1);
#pragma omp parallel num_threads(2) 
    {
      if(omp_get_thread_num() == 0)
	learner->Train(strlen(outfile) ? outfile : NULL, true);
      else 
	retval = RunServer(port);
    }
  } else
    retval = RunServer(port);
  return retval;
}



int StructuredLearnerRpc::RunServer(int port) {
  AddMethods();
  server->PrintUsage();
  server->RunServer();

  return 0;
}



// Add remote methods for online learning, classification, and active labeling, and document their parameters/return values
void StructuredLearnerRpc::AddMethods() {
  Json::Value new_session_parameters, new_session_returns;
  new_session_parameters["x"] = "String encoding of an example in the format of StructuredData::load()";
  new_session_returns["session_id"] = "A string encoding of the new session id.  The client should pass this as a parameter to all future accesses to x";
  server->RegisterMethod(new JsonRpcMethod<StructuredLearnerRpc>(this, &StructuredLearnerRpc::NewSession, "new_session", "Start a new session, where a session caches data for processing an example over a series of operations over the network", new_session_parameters, new_session_returns));

  Json::Value classify_example_parameters, classify_example_returns;
  classify_example_parameters["session_id"] = "Optional session id returned by new_session().  If this doesn't exist, a new session is created.";
  classify_example_parameters["x"] = "Optional string encoding of an example in the format of StructuredData::load().  If this was already passed to a previous request on this session, it is unnecessary.";
  classify_example_parameters["partial_label"] = "Optional string encoding of a partial label in the format of StructuredLabel::load().  The returned label y must be consistent with the partial label.";
  classify_example_returns["y"] = "A string encoding of the predicted label y in the format of StructuredLabel::load()";
  classify_example_returns["score"] = "The score of the predicted label y";
  classify_example_returns["session_id"] = "A string encoding of the session id.  The client should pass this as a parameter to all future accesses to x";
  server->RegisterMethod(new JsonRpcMethod<StructuredLearnerRpc>(this, &StructuredLearnerRpc::ClassifyExample, "classify_example", "Classify an example: given an example x predict the label y with highest score", classify_example_parameters, classify_example_returns));


  if(train) {
    Json::Value add_example_parameters, add_example_returns;
    add_example_parameters["session_id"] = "Optional session id returned by new_session().";
    add_example_parameters["y"] = "String encoding of the ground truth example label format of StructuredLabel::load()";
    add_example_parameters["x"] = "Optional string encoding of an example in the format of StructuredData::load().  If this was already passed to a previous request on this session, it is unnecessary.";
    add_example_returns["index"] = "An integer index into the training set list of examples";
    add_example_returns["session_id"] = "A string encoding of the session id.  The client should pass this as a parameter to all future accesses to x";
    server->RegisterMethod(new JsonRpcMethod<StructuredLearnerRpc>(this, &StructuredLearnerRpc::AddNewExample, "add_example", "Add a new training example to the structured learner, which is currently training in online fashion", add_example_parameters, add_example_returns));

    //server->RegisterMethod(new JsonRpcMethod<StructuredLearnerRpc>(*this, &StructuredLearnerRpc::GetStatistics, "plot_stats", get_statistics));

    Json::Value save_parameters;
    save_parameters["filename"] = "The filename defining where to save the file";
    save_parameters["saveFull"] = "If true, saves sufficient information to resume online learning again.  Otherwise, only saves enough info to classify examples.";
    server->RegisterMethod(new JsonRpcMethod<StructuredLearnerRpc>(this, &StructuredLearnerRpc::SaveCurrentModel, "save", "Save all data associated with the online learner to disk", save_parameters));

    server->RegisterMethod(new JsonRpcMethod<StructuredLearnerRpc>(this, &StructuredLearnerRpc::Shutdown, "shutdown", "Signal the server to shutdown"));


    Json::Value set_parameter_parameters;
    set_parameter_parameters["C"] = "Optional regularization parameter C (C=1/lambda).";
    set_parameter_parameters["lambda"] = "Optional regularization parameter lambda (C=1/lambda).";
    set_parameter_parameters["featureScale"] = "Optional parameter that scales the feature space Psi(x,y).";
    server->RegisterMethod(new JsonRpcMethod<StructuredLearnerRpc>(this, &StructuredLearnerRpc::SetParameter, "set_parameter", "Set a learning parameter, which can be one of C, lambda, or featureScale", set_parameter_parameters));
  }
   
  Json::Value evaluate_testset_parameters;
  evaluate_testset_parameters["testset"] = "Filename of the testset, in the format of StructuredSVM::LoadDataset()";
  evaluate_testset_parameters["predictions"] = "Filename of where to store predictions, where each line of the file corresponds to one test example in the format <y_predicted> <y_groundtruth> <loss> <score_predicted> <score_groundtruth>, and <y_predicted> and <y_groundtruth> are in the format of StructuredLabel::load()";
  server->RegisterMethod(new JsonRpcMethod<StructuredLearnerRpc>(this, &StructuredLearnerRpc::EvaluateTestset, "evaluate_testset", "Evaluate performance on a testset", evaluate_testset_parameters));
}
