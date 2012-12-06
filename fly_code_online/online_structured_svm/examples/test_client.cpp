#include <cstdio>
#include <cstdlib>

/**
 * @example test_client.cpp
 *
 * Simple client program for simulated online interactive training.  It is used only to give an example of
 * how the network protocol works and is otherwise not useful.  The server can be started
 * initially with no training examples.  This client will connect into the server over
 * a network socket, and incrementally add new training examples.  Before adding a new
 * example, the client asks the server to predict its label using the current learned
 * model, and we record the average loss
 *
 * Example usage:
 *   - Start the server: 
\htmlonly <div style="padding: 0.5em 1em; border-top: 1px solid #ddd; border-bottom: 1px solid #ddd; background-color: #eaeafa;">
$ examples/bin/release_static/structured_svm_multiclass.out -p classes.txt -P 8096
</div> \endhtmlonly
 *   - Start the client: 
\htmlonly <div style="padding: 0.5em 1em; border-top: 1px solid #ddd; border-bottom: 1px solid #ddd; background-color: #eaeafa;">
$ examples/bin/release_static/test_client.out train.txt 8096
</div> \endhtmlonly
 *
 */

#include "structured_svm_multiclass.h"
#include "jsonrpc.h"

#define IP "127.0.0.1"   /**< localhost */
#define TARGET_TIME 100  /**< Make the simulated training/labeling session last about 100 seconds */


int main(int argc, char** argv)
{
  /* avoid compilation warnings */
  argc = argc;
  argv = argv;

  if(argc < 3) {
    fprintf(stderr, "USAGE: ./test_client.out <train_file> <port>\n");
    return -1;
  }

  char trainfile[1000];

  int id = 1;
  MulticlassStructuredSVM params;
  strcpy(trainfile, argc > 1 ? argv[1] : "data/train.dat");
  int port = atoi(argv[2]);
  
  if(argc > 3)
    assert(params.Load(argv[3]));


  StructuredDataset *trainset = params.LoadDataset(trainfile);
  double sumLoss = 0;

  for(int i = 0; i < trainset->num_examples; i++) {

    /* build JSON-RPC query */
    Json::Value query1;
    query1["jsonrpc"] = "2.0";
    query1["id"] = id++;
    query1["method"] = "classify_example";
    query1["x"] = trainset->examples[i]->x->save(&params);

    Json::Value response = JsonRpcClientRequest(IP, port, query1);
    if(!response.isMember("session_id")) {
      std::cerr << "Didn't get session_id response from server!" << std::endl;
      exit(EXIT_FAILURE);
    }
    if(!response.isMember("y")) {
      std::cerr << "Didn't get label response from server!" << std::endl;
      exit(EXIT_FAILURE);
    }

    /*******************************************************************
     * Here, one can present response["label"] to the user, and
     * make a sequence of subsequent requests of the form:
     *   query1["method"] = "classify_example";
     *   query1["partial_label"] = <user_specified_partial_labeling>
     *   query1["session_id"] = response["session_id"];
     *******************************************************************/

    StructuredLabel *y = params.NewStructuredLabel(trainset->examples[i]->x);
    assert(y->load(response["y"], &params));
    Json::FastWriter writer;
    double l = params.Loss(trainset->examples[i]->y, y);
    std::cout << "  " << i << ": predicted " << writer.write(response["y"]) << 
      " when true label is " << writer.write(trainset->examples[i]->y->save(&params)) << " (loss=" << l << ")\n";
    sumLoss += l;
    delete y;

    Json::Value query2;
    query2["jsonrpc"] = "2.0";
    query2["id"] = id++;
    query2["method"] = "add_example";
    query2["session_id"] = response["session_id"];
    query2["y"] = trainset->examples[i]->y->save(&params);
    response = JsonRpcClientRequest(IP, port, query2);

    usleep(TARGET_TIME * 1000000 / trainset->num_examples);
  }

  std::cout << "Average loss was " << (sumLoss / trainset->num_examples) << "\n";


  return EXIT_SUCCESS;
}

