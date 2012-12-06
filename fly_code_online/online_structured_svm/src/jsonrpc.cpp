#include "jsonrpc.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifdef WIN32
#include <Winsock2.h>
#include <WS2tcpip.h>
#include <io.h>
#else
#include <unistd.h>
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
  SOCKET sockfd = socket(AF_INET, SOCK_STREAM, 0);
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
    SOCKET newsockfd = accept(sockfd, (struct sockaddr *)&cli_addr, &clilen);

    if(shutdown) break;

    if (newsockfd < 0) { 
      fprintf(stderr, "ERROR accepting socket connection\n");
      return;
    }

#ifdef WIN32
    int n = recv(newsockfd, buf, 1999999, 0);
#else
    int n = read(newsockfd, buf, 1999999);
#endif
    buf[n] = 0;

    ProcessRequest(newsockfd, buf);
  }
#ifdef WIN32
  closesocket(sockfd);
  WSACleanup( );
#else
  close(sockfd);
#endif
}
  
void JsonRpcServer::ProcessRequest(SOCKET  newsockfd, char *buf) {
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

#ifdef WIN32
  closesocket(newsockfd);
#else
  close(newsockfd);
#endif
}

void JsonRpcServer::ReturnJsonResponse(SOCKET  newsockfd, const Json::Value &response) {
  Json::FastWriter writer;
  char *str = (char*)malloc(500000);
  strcpy(str, writer.write(response).c_str());
#ifdef WIN32
  send(newsockfd, str, strlen(str), 0);
#else
  write(newsockfd, str, strlen(str));
#endif
  free(str);
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
  SOCKET sock = socket (AF_INET, SOCK_STREAM, 0);
  if (sock == -1) {
    fprintf(stderr, "Failed to create socket\n");
    return Json::Value::null;
  }
  sin.sin_family = AF_INET;
  sin.sin_port = htons( (unsigned short)port);

  struct hostent *host_addr = gethostbyname(hostname);
  if(host_addr==NULL) {
#ifdef WIN32
    closesocket(sock);
#else
    close(sock);
#endif
    fprintf(stderr, "Failed to find hostname %s\n", hostname);
    return Json::Value::null;
  }
  sin.sin_addr.s_addr = *((int*)*host_addr->h_addr_list) ;

  if( connect (sock,(const struct sockaddr *)&sin, sizeof(sockaddr_in) ) == -1 ) {
#ifdef WIN32
    closesocket(sock);
#else
    close(sock);
#endif
    fprintf(stderr, "Failed to find connect to socket\n");
    return Json::Value::null;
  }

  Json::FastWriter writer;
  char *str = (char*)malloc(500000);
  strcpy(str, writer.write(req).c_str());
  send(sock,str,strlen(str),0);
  int n = recv(sock,str,100000,0);
  str[n] = '\0';
  //fprintf(stderr, "%s\n", str);
  

#ifdef WIN32
  closesocket(sock);
  WSACleanup( );
#else
  close(sock);
#endif

  Json::Reader reader;
  Json::Value ret;
  if(!reader.parse(str, ret))
	  ret = Json::Value::null;
  free(str);
  return ret;
}
