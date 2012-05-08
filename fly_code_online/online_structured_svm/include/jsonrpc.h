#ifndef __JSONRPC_H
#define __JSONRPC_H

#include "json/json.h"
#include "util.h"

class JsonRpcMethodAbstract {
public:
  virtual bool Invoke(const Json::Value& msg, Json::Value& response) = 0;
  virtual Json::Value Usage() = 0;
  virtual const char *Name() = 0;
  virtual ~JsonRpcMethodAbstract() {}
};

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

  bool Invoke(const Json::Value& msg, Json::Value& response){
    return (server->*callback)(msg, response);
  }

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
  
  const char *Name() { return name; }
};

class JsonRpcServer {
  JsonRpcMethodAbstract **functions;
  int num_functions;
  int port;
  bool shutdown;

 public: 
  JsonRpcServer(int port);
  ~JsonRpcServer();
  void RegisterMethod(JsonRpcMethodAbstract *func);
  void UnregisterMethod(const char *name);
  void RunServer();
  void ProcessRequest(int newsockfd, char *buf);
  void ReturnJsonResponse(int newsockfd, const Json::Value &response);
  JsonRpcMethodAbstract *FindMethod(const char *name);
  int FindMethodInd(const char *name);
  Json::Value Usage();
  void PrintUsage();
  void Shutdown();
};

Json::Value JsonRpcClientRequest(const char *hostname, int port, Json::Value &req);

#endif
