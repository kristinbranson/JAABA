#ifndef __JSONRPC_H
#define __JSONRPC_H

#include "json/json.h"
#include "util.h"

/**
 * @file jsonrpc.h
 * @brief Implements a JSON-RPC server (interface for hosting functions that can be invoked over the network using a JSON protocol)
 */

#ifdef WIN32
#include <Winsock2.h>
#include <WS2tcpip.h>
#include <io.h>
#else
#include <unistd.h>
#include <netdb.h>
#endif

#ifndef WIN32
/// @cond
#define SOCKET int
/// @endcond
#endif



class StructuredLearnerRpc;

/// @cond
class JsonRpcMethodAbstract {
public:
  virtual bool Invoke(const Json::Value& msg, Json::Value& response) = 0;
  virtual Json::Value Usage() = 0;
  virtual const char *Name() = 0;
  virtual ~JsonRpcMethodAbstract() {}
};
/// @endcond

/**
 * @class JsonRpcMethod
 * @brief Container class for a function that is callable as a remote procedure call over the network
 */
template<class T> 
class JsonRpcMethod : public JsonRpcMethodAbstract {
  T *server;
  char *name;
  char *description;
  Json::Value parameters;
  Json::Value return_values;
  bool (T::*callback)(const Json::Value& msg, Json::Value& response);

  friend class JsonRpcServer;

public:
  /**
   * @brief Create a new method
   * @param server Pointer to an object of some class that inherits from StructuredLearnerRpc and defines a bunch of member functions
   * that are callable as remote procedure calls
   * @param callback a pointer to a member function of the same class as server
   * @param name The name of the method (this is the name that people will use to invoke this procedure over the network)
   * @param description If non-null, a description of the method that will be displayed when printing the USAGE options for this function
   * @param parameters If non-null, a JSON encoding of an array of possible parameters that will be displayed when printing the USAGE options for this function
   * @param return_values If non-null, a JSON encoding of an array of possible return values that will be displayed when printing the USAGE options for this function
   */
  JsonRpcMethod(T *server,  bool (T::*callback)(const Json::Value& msg, Json::Value& response),
		const char *name, const char *description = NULL, Json::Value parameters  = Json::Value::null,
		Json::Value return_values  = Json::Value::null) {
    this->server = server;
    this->callback = callback;
    this->name = StringCopy(name);
    this->description = description ? StringCopy(description) : NULL;
    this->parameters = parameters;
    this->return_values = return_values;
  }

  ~JsonRpcMethod() {
    if(name) free(name);
    if(description) free(description);
  }

  /**
   * @brief Invoke the function for this method
   * @param msg A JSON object storing an array of parameters
   * @param response A JSON object into which an array of return values is written
   */
  bool Invoke(const Json::Value& msg, Json::Value& response){
    return (server->*callback)(msg, response);
  }

  /**
   * @brief Print the usage options giving a description of possible parameters and return values for this function
   */
  Json::Value Usage() {
    if(!description) return Json::Value::null;
    else {
      Json::Value retval;
      if(description) retval["description"] = description;
      if(parameters != Json::Value::null) retval["parameters"] = parameters;
      if(return_values != Json::Value::null) retval["return_values"] = return_values;
      return retval;
    }
  }
  
  /**
   * @brief Get the name of this method
   */
  const char *Name() { return name; }
};

/**
 * @class JsonRpcServer
 * @brief An abstract class that handles network communication for implementing remote procedure calls via a JSON-RPC protocol.  
 */
class JsonRpcServer {
  friend class StructuredLearnerRpc;
  JsonRpcMethodAbstract **functions;
  int num_functions;
  int port;
  bool shutdown;


  JsonRpcServer(int port);
  ~JsonRpcServer();
  void RunServer();
  void ProcessRequest(SOCKET  newsockfd, char *buf);
  void ReturnJsonResponse(SOCKET  newsockfd, const Json::Value &response);
  JsonRpcMethodAbstract *FindMethod(const char *name);
  int FindMethodInd(const char *name);
  Json::Value Usage();
  void PrintUsage();
  void Shutdown();

 public:

  /**
   * @brief Register a new method to be available as a remote procedure call
   * @param func A JsonRpcMethod of the function to add
   */
  void RegisterMethod(JsonRpcMethodAbstract *func);

  /**
   * @brief Remove a previously registered method
   * @param name The name of a JsonRpcMethod to remove
   */
  void UnregisterMethod(const char *name);
};

/**
 * @brief Connect into a server over the network to invoke a function using a remote procedure call
 * @param hostname The ip or host name of the server
 * @param port The network port of the server
 * @param req A JSON object encoding the appropriate parameters and headers for a JSON-RPC client request
 */
Json::Value JsonRpcClientRequest(const char *hostname, int port, Json::Value &req);

#endif
