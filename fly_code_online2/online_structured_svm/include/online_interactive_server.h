#ifndef __SERVER_H
#define __SERVER_H

#include "structured_svm.h"
#include "jsonrpc.h"

#include <json/json.h>
#include <cstdio>
#include <cstdlib>
#include <csignal>

#include <omp.h>

#define MAX_SESSIONS 20  /**< The maximum number of network sessions that are kept cached in memory */

/**
 * @file online_interactive_server.h
 * @brief Implements a network server for training a structured SVM while allowing new examples to be added in online fashion or labeled interactively
 */




/**
 * @struct Session
 * @brief Caches memory for a network session.  
 *
 * For example, if a person is interactively
 *   classifying an example x, it may be convenient to keep x cached in memory as the
 *   client sends repeated commands over the network
 */
typedef struct {
  StructuredExample *example; /**< A memory cache of the example */
  StructuredLabel *partial_label;  /**< The last partial label specified by the user */

  char *id;   /**< a unique id for this session */
  unsigned long timestamp;  /**< a timestep when this session was last accessed */
  omp_lock_t lock;  /**< synchronizes access to this session's memory */
} Session;


/**
 * @class StructuredLearnerRpc
 * @brief A class used to relay commands from a client over the network to an online structured
 *   learner in the process of training.  This is useful for dynamically adding new examples or
 *   interactively classifying examples while learning is in progress
 * 
 * In the documentation below, we describe the network protocol for each type of client request,
 * where "Client: client_str" means the client sends a string client_str to the server over a network
 * socket connection, and "Server: server_str" means the server sends a response string server_str back.
 * Each type of client request results in invoking one of the member functions of this class
 * with the applicable arguments.  To test this, you can do something like this from a unix prompt:
\htmlonly <div style="padding: 0.5em 1em; border: 1px solid #ccf; background-color: #eaeafa;">
   $ telnet localhost 8086
   <br>&nbsp;&nbsp;client_str
</div> \endhtmlonly 
 *
 */
class StructuredLearnerRpc
{
protected:
/// @cond
  StructuredSVM *learner;
  unsigned long timestamp;
  JsonRpcServer *server;
  
  Session sessions[MAX_SESSIONS];  /**< a num_sessions array of pointers to sessions */
  int num_sessions;

  char session_dir[1000];
  omp_lock_t lock;

  char trainfile[1000], infile[1000], initial_sample_set[1000];
  bool train;
  bool randomize;
/// @endcond


public:

  /**
   * @brief Create a new StructuredLearnerRpc, where l defines the class used for structured learning
   */
  StructuredLearnerRpc(StructuredSVM *l);
  ~StructuredLearnerRpc();

  /**
   * @brief Begin running a server
   */
  int main(int argc, const char **argv) ;

protected:
  /**
   * @brief Add remote procedure methods
   *
   */
  virtual void AddMethods();

  /**
   * @brief Parse command line arguments
   *
   */
  virtual void parse_command_line_arguments(int argc, const char **argv);

  /**
   * @brief Add a new example to the training set
   * @param root A JSON object storing an array of parameters
   * @param response A JSON object into which an array of return values is written
   *
   * The network protocol for adding a new example (x,y) where json_encoding_of_x and json_encoding_of_y are in
   * the format of StructuredData::write() and StructuredLabel::write() and ind is the index of the training example added:
\htmlonly <div style="padding: 0.5em 1em; border: 1px solid #aaf; background-color: #dadafa;">
<font color="blue">Client:</font> {"jsonrpc":"2.0","id":1,"method":"add_example","x":json_encoding_of_x,"y":json_encoding_of_y}
</div> \endhtmlonly
\htmlonly <div style="padding: 0.5em 1em; border: 1px solid #faa; background-color: #fadada;">
<font color="red">Server:</font> {"id":1,"session_id":"session_id","index":ind}
</div> \endhtmlonly
   *
   */
  virtual bool AddNewExample(const Json::Value& root, Json::Value& response);

  /**
   * @brief Change the label of an existing training example
   * @param root A JSON object storing an array of parameters
   * @param response A JSON object into which an array of return values is written
   *
   * The network protocol for changing the label example, where json_encoding_of_y is in
   * the format of StructuredLabel::write() and example_index is the index of the training 
   * example to edit (as returned by AddNewExample()):
\htmlonly <div style="padding: 0.5em 1em; border: 1px solid #aaf; background-color: #dadafa;">
<font color="blue">Client:</font> {"jsonrpc":"2.0","id":1,"method":"relabel_example","index":example_index,"y":json_encoding_of_y}
</div> \endhtmlonly
\htmlonly <div style="padding: 0.5em 1em; border: 1px solid #faa; background-color: #fadada;">
<font color="red">Server:</font> {"id":1,"session_id":"session_id","index":ind}
</div> \endhtmlonly
   *
   */
  virtual bool RelabelExample(const Json::Value& root, Json::Value& response);

  /**
   * @brief Classify an example x, which is specified by the client.  The server returns a predicted label y.  Takes an optional partial label
   * @param root A JSON object storing an array of parameters
   * @param response A JSON object into which an array of return values is written
   *
   *  The network protocol for a regular classification is (where json_encoding_of_x is in the format
   *  of StructuredData::write() and is specified by the client, and json_encoding_of_y is in the
   *  format of StructuredLabel::write() and is generated by the server:
\htmlonly <div style="padding: 0.5em 1em; border: 1px solid #aaf; background-color: #dadafa;">
<font color="blue">Client:</font> {"jsonrpc":"2.0","id":1,"method":"classify_example","x":"json_encoding_of_x"}
</div> \endhtmlonly
\htmlonly <div style="padding: 0.5em 1em; border: 1px solid #faa; background-color: #fadada;">
<font color="red">Server:</font> {"id":1,"session_id":"session_id","y":"json_encoding_of_y"}
</div> <br> \endhtmlonly
   *
   * One can optionally do a sequence of classification requests on the same example, for interactive classification. The following example
   *  does 3 steps of interactive labeling then adds the example to the training set:
   * - Client sends a new example json_encoding_of_x in the format of StructuredData::write(), and the server sends back the initial max likelihood solution json_encoding_of_y1 in the format of StructuredLabel::write()
\htmlonly <div style="padding: 0.5em 1em; border: 1px solid #aaf; background-color: #dadafa;">
<font color="blue">Client:</font> {"jsonrpc":"2.0","id":1,"method":"classify_example","x":"json_encoding_of_x"}
</div> \endhtmlonly
\htmlonly <div style="padding: 0.5em 1em; border: 1px solid #faa; background-color: #fadada;">
<font color="red">Server:</font> {"id":1,"session_id":"session_id","y":"json_encoding_of_y1"}
</div> <br> \endhtmlonly
 *
 * - User provides a partial labeling partial_label1 in the format of StructuredLabel::write() that corrects one of the variables in y, and the server sends back the maximum likelihood solution json_encoding_of_y2 that is consistent with the partial labeling in the format of StructuredLabel::write()
\htmlonly <div style="padding: 0.5em 1em; border: 1px solid #aaf; background-color: #dadafa;">
<font color="blue">Client:</font> {"jsonrpc":"2.0","id":1,"method":"classify_example","session_id":"session_id","partial_label":"partial_label1"}
</div> \endhtmlonly
\htmlonly <div style="padding: 0.5em 1em; border: 1px solid #faa; background-color: #fadada;">
<font color="red">Server:</font> {"id":1,"session_id":"session_id","y":"json_encoding_of_y2"}
</div> <br> \endhtmlonly
 *
 * - User provides a partial labeling partial_label2 in the format of StructuredLabel::write() that corrects another variable in y, and the server sends back the maximum likelihood solution json_encoding_of_y3 that is consistent with the partial labeling in the format of StructuredLabel::write()
\htmlonly <div style="padding: 0.5em 1em; border: 1px solid #aaf; background-color: #dadafa;">
<font color="blue">Client:</font> {"jsonrpc":"2.0","id":1,"method":"classify_example","session_id":"session_id","partial_label":"partial_label2"}
</div> \endhtmlonly
\htmlonly <div style="padding: 0.5em 1em; border: 1px solid #faa; background-color: #fadada;">
<font color="red">Server:</font> {"id":1,"session_id":"session_id","y":"json_encoding_of_y3"}
</div> <br> \endhtmlonly
 *
 * - When the user is satisfied with the solution json_encoding_of_y3, the client asks the server to add a new verified training example, and the server sends back its index
\htmlonly <div style="padding: 0.5em 1em; border: 1px solid #aaf; background-color: #dadafa;">
<font color="blue">Client:</font> {"jsonrpc":"2.0","id":1,"method":"add_example","session_id":"session_id","y":"json_encoding_of_y3"}
</div> \endhtmlonly
\htmlonly <div style="padding: 0.5em 1em; border: 1px solid #faa; background-color: #fadada;">
<font color="red">Server:</font> {"id":1,"session_id":"session_id","index":ind}
</div> <br> \endhtmlonly
   *
   */
  virtual bool ClassifyExample(const Json::Value& root, Json::Value& response);

  /**
   * @brief Save the current model for a training job in progress
   * @param root A JSON object storing an array of parameters
   * @param response A JSON object into which an array of return values is written
   *
   * The network protocol is:
\htmlonly <div style="padding: 0.5em 1em; border: 1px solid #aaf; background-color: #dadafa;">
<font color="blue">Client:</font> {"jsonrpc":"2.0","id":1,"method":"save","filename":"filename_to_save_to"}
</div> \endhtmlonly
\htmlonly <div style="padding: 0.5em 1em; border: 1px solid #faa; background-color: #fadada;">
<font color="red">Server:</font> {"id":1}
</div> <br> \endhtmlonly
   */
  virtual bool SaveCurrentModel(const Json::Value& root, Json::Value& response);

  /**
   * @brief Signal the current training algorithm to exit
   * @param root A JSON object storing an array of parameters
   * @param response A JSON object into which an array of return values is written
   *
   * The network protocol is:
\htmlonly <div style="padding: 0.5em 1em; border: 1px solid #aaf; background-color: #dadafa;">
<font color="blue">Client:</font> {"jsonrpc":"2.0","id":1,"method":"shutdown"}
</div> \endhtmlonly
\htmlonly <div style="padding: 0.5em 1em; border: 1px solid #faa; background-color: #fadada;">
<font color="red">Server:</font> {"id":1}
</div> <br> \endhtmlonly
   */
  virtual bool Shutdown(const Json::Value& root, Json::Value& response);

  /**
   * @brief Change a learning parameter C, lambda, or feature_scale
   * @param root A JSON object storing an array of parameters
   * @param response A JSON object into which an array of return values is written
   *
   * The network protocol is:
\htmlonly <div style="padding: 0.5em 1em; border: 1px solid #aaf; background-color: #dadafa;">
<font color="blue">Client:</font> {"jsonrpc":"2.0","id":1,"method":"set_parameter","C":c}
</div> \endhtmlonly
\htmlonly <div style="padding: 0.5em 1em; border: 1px solid #faa; background-color: #fadada;">
<font color="red">Server:</font> {"id":1}
</div> <br> \endhtmlonly
   */
  virtual bool SetParameter(const Json::Value& root, Json::Value& response);

  /**
   * @brief Evaluate a testset using the current model parameters
   * @param root A JSON object storing an array of parameters
   * @param response A JSON object into which an array of return values is written
   *
   * The network protocol is:
\htmlonly <div style="padding: 0.5em 1em; border: 1px solid #aaf; background-color: #dadafa;">
<font color="blue">Client:</font> {"jsonrpc":"2.0","id":1,"method":"evaluate_testset","testset":"testset_filename"}
</div> \endhtmlonly
\htmlonly <div style="padding: 0.5em 1em; border: 1px solid #faa; background-color: #fadada;">
<font color="red">Server:</font> {"id":1,"ave_loss":loss}
</div> <br> \endhtmlonly
   */
  virtual bool EvaluateTestset(const Json::Value& root, Json::Value& response);

  /**
   * @brief Plot a decomposition of the test error
   * @param root A JSON object storing an array of parameters
   * @param response A JSON object into which an array of return values is written
   *
   * The network protocol is:
\htmlonly <div style="padding: 0.5em 1em; border: 1px solid #aaf; background-color: #dadafa;">
<font color="blue">Client:</font> {"jsonrpc":"2.0","id":1,"method":"plot_stats","plot_name":"plot_by_time","plot_by":"time"}
</div> \endhtmlonly
\htmlonly <div style="padding: 0.5em 1em; border: 1px solid #faa; background-color: #fadada;">
<font color="red">Server:</font> {"id":1}
</div> <br> \endhtmlonly
   * where one can use "time", "example", or "iteration" for the "plot_by" parameter 
   */
  virtual bool GetStatistics(const Json::Value& root, Json::Value& response);

  /**
   * @brief Create a new network session (which caches in memory calculations pertaining to a particular example x
   * @param root A JSON object storing an array of parameters
   * @param response A JSON object into which an array of return values is written
   *
   * The network protocol is, where json_encoding_of_x is in the format of StructuredData:write():
\htmlonly <div style="padding: 0.5em 1em; border: 1px solid #aaf; background-color: #dadafa;">
<font color="blue">Client:</font> {"jsonrpc":"2.0","id":1,"method":"new_session","x":"json_encoding_of_x"}
</div> \endhtmlonly
\htmlonly <div style="padding: 0.5em 1em; border: 1px solid #faa; background-color: #fadada;">
<font color="red">Server:</font> {"id":1,"session_id":"session_id"}
</div> <br> \endhtmlonly
   */
  virtual bool NewSession(const Json::Value& root, Json::Value& response);

  virtual bool Visualize(const Json::Value& root, Json::Value& response);

  virtual bool CheckAlive(const Json::Value& root, Json::Value& response) {
    response["message"] = "hello";
    return true;
  }

/// @cond
 protected:
  int FindOrCreateSession(const Json::Value& root, Json::Value& response);

  bool InitializeSession(const Json::Value& root, Json::Value& response);
  
  int FindSession(const char *sess_id, bool lockSession=false);
  int FindSession(const Json::Value& root, Json::Value& response);
  void EvictSession(int ind=-1);
  void Lock() { omp_set_lock(&lock); }
  void Unlock() { omp_unset_lock(&lock); }
  void LockSession(int i) { omp_set_lock(&sessions[i].lock); }
  void UnlockSession(int i) { omp_unset_lock(&sessions[i].lock); }

  int RunServer(int port);

 protected:
  char paramfile[1000], testfile[1000], outfile[1000], predictionsfile[1000], plotfile[1000];
  int port;
  bool runServer;
/// @endcond
};

/// @cond
#define JSON_ERROR(errStr, sess_ind) \
  fprintf(stderr, "%s\n", errStr); \
  response["error"] = std::string(errStr); \
  if(sess_ind >= 0) UnlockSession(sess_ind); \
  return false;    
/// @endcond


#endif
