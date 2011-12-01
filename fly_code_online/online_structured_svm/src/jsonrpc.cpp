#include "jsonrpc.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#ifdef WIN32
#include <Winsock2.h>
#else
#include <netdb.h> 
#endif


JsonRpcServer::JsonRpcServer(int port) {
  functions = NULL;
  num_functions = 0;
  this->port = port;
}

JsonRpcServer::~JsonRpcServer() {
  for(int i = 0; i < num_functions; i++)
    delete functions[i];
  free(functions);
}

JsonRpcMethodAbstract *JsonRpcServer::FindMethod(const char *name) {
  int ind = FindMethodInd(name);
  return ind >= 0 ? functions[ind] : NULL;
}

int JsonRpcServer::FindMethodInd(const char *name) {
  for(int i = 0; i < num_functions; i++)
    if(!strcmp(functions[i]->Name(), name))
      return i;
  return -1;
}

void JsonRpcServer::RegisterMethod(JsonRpcMethodAbstract *func) {
  if(!FindMethod(func->Name())) {
    functions = (JsonRpcMethodAbstract**)realloc(functions, sizeof(JsonRpcMethodAbstract*)*(num_functions+1));
    functions[num_functions++] = func;
  }
}

void JsonRpcServer::UnregisterMethod(const char *name) {
  int ind = FindMethodInd(name);
  if(ind >= 0) {
    for(int i = ind; i < num_functions-1; i++)
      functions[i] = functions[i+1];
    delete functions[ind];
    num_functions--;
  }
}

void JsonRpcServer::Shutdown() {
  shutdown = true;
  Json::Value req;
  req["method"] = "terminate";
  JsonRpcClientRequest("127.0.0.1", port, req);
}

void JsonRpcServer::RunServer() {
#ifdef WIN32
  {
    WSADATA WsaData;
    WSAStartup (0x0101, &WsaData);
  }
#endif

  shutdown = false;

  struct sockaddr_in serv_addr, cli_addr;
  int sockfd = socket(AF_INET, SOCK_STREAM, 0);
  if (sockfd < 0) {
    fprintf(stderr, "ERROR opening socket\n");
    return;
  }
  
  memset((char *) &serv_addr, 0, sizeof(serv_addr));
  serv_addr.sin_family = AF_INET;
  serv_addr.sin_addr.s_addr = INADDR_ANY;
  serv_addr.sin_port = htons(port);
  if (bind(sockfd, (struct sockaddr *) &serv_addr, sizeof(serv_addr)) < 0) {
    fprintf(stderr, "ERROR binding socket\n");
    return;
  }


  fprintf(stderr, "Server listening on port %d\n", port);

  listen(sockfd,5);
  socklen_t clilen = sizeof(cli_addr);

  char *buf = new char[2000000];
  while(!shutdown) {
    int newsockfd = accept(sockfd, (struct sockaddr *)&cli_addr, &clilen);

    if(shutdown) break;

    if (newsockfd < 0) { 
      fprintf(stderr, "ERROR accepting socket connection\n");
      return;
    }

    int n = read(newsockfd, buf, 1999999);
    buf[n] = 0;

    ProcessRequest(newsockfd, buf);
  }
  close(sockfd);

 #ifdef WIN32
  WSACleanup( );
 #endif

}
  
void JsonRpcServer::ProcessRequest(int newsockfd, char *buf) {
  Json::Value root;
  Json::Value response;
  Json::Reader reader;

  response["jsonrpc"] = "2.0";
  if(root.isMember("id")) response["id"] = root["id"];
  
  if(!reader.parse(buf, root)) {
    response["error"] = "Error parsing request";  ReturnJsonResponse(newsockfd, response);
  } else {
    if(!root.isMember("jsonrpc") || root.get("jsonrpc", "").asString() != "2.0") {
      response["error"] = "Error parsing request, parameter \"jsonrpc\" must be \"2.0\""; 
      ReturnJsonResponse(newsockfd, response);
    } else if(!root.isMember("method")) {
      response["error"] = "Error parsing request, parameter \"method\" not found"; 
      ReturnJsonResponse(newsockfd, response);
    } else {
      char name[1000];
      strcpy(name, root.get("method", "").asString().c_str());
      JsonRpcMethodAbstract *func = FindMethod(name);
      if(!func) {
	char str[1000]; sprintf(str, "Error parsing request, \"method\" %s not found", name);
	response["error"] = str; 
	ReturnJsonResponse(newsockfd, response);
      } else {
	// TODO: Invoke this as a thread
	func->Invoke(root, response);
	ReturnJsonResponse(newsockfd, response);
      }
    }
  }
  close(newsockfd);
}

void JsonRpcServer::ReturnJsonResponse(int newsockfd, const Json::Value &response) {
  Json::FastWriter writer;
  char str[500000];
  strcpy(str, writer.write(response).c_str());
  write(newsockfd, str, strlen(str));
}
  

Json::Value JsonRpcServer::Usage() {
  Json::Value retval;
  for(int i = 0; i < num_functions; i++) {
    Json::Value u = functions[i]->Usage();
    if(u != Json::Value::null) 
      retval[functions[i]->Name()] = u;
  }
  return retval;
}

void JsonRpcServer::PrintUsage() {
  Json::Value usage = Usage();
  Json::StyledWriter writer;
  fprintf(stderr, "%s\n", writer.write(usage).c_str());
}



Json::Value JsonRpcClientRequest(const char *hostname, int port, Json::Value &req) {
#ifdef WIN32
  {
    WSADATA WsaData;
    WSAStartup (0x0101, &WsaData);
  }
#endif
  
  req["jsonrpc"] = "2.0";
  req["id"] = "0";


  sockaddr_in sin;
  int sock = socket (AF_INET, SOCK_STREAM, 0);
  if (sock == -1) {
    fprintf(stderr, "Failed to create socket\n");
    return Json::Value::null;
  }
  sin.sin_family = AF_INET;
  sin.sin_port = htons( (unsigned short)port);

  struct hostent *host_addr = gethostbyname(hostname);
  if(host_addr==NULL) {
    close(sock);
    fprintf(stderr, "Failed to find hostname %s\n", hostname);
    return Json::Value::null;
  }
  sin.sin_addr.s_addr = *((int*)*host_addr->h_addr_list) ;

  if( connect (sock,(const struct sockaddr *)&sin, sizeof(sockaddr_in) ) == -1 ) {
    close(sock);
    fprintf(stderr, "Failed to find connect to socket\n");
    return Json::Value::null;
  }

  Json::FastWriter writer;
  char str[100000];
  strcpy(str, writer.write(req).c_str());
  send(sock,str,strlen(str),0);
  int n = recv(sock,str,100000,0);
  str[n] = '\0';
  //fprintf(stderr, "%s\n", str);
  
  close(sock);

#ifdef WIN32
  WSACleanup( );
#endif

  Json::Reader reader;
  Json::Value ret;
  return reader.parse(str, ret) ? ret : Json::Value::null;
}
