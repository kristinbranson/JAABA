#include <cstdio>
#include <cstdlib>
#include <assert.h>
#include <stdio.h>

#include <jsonrpc/jsonrpc.h>


// USAGE: ./client commands.json
//        Executes the commands in the file commands.json
//
//        or
//        
//        ./client 
//        Reads json commands from stdin

void SendJSONCommand(std::string queryStr);


int main(int argc, char** argv)
{
  if(!networking::init()) {
      std::cerr << "Networking initialization failed" << std::endl;
      exit(EXIT_FAILURE);
  }

  const char *commandFile = argc <= 1 ? NULL : argv[1];
	
  FILE *fin = commandFile ? fopen(argv[1], "r") : stdin;
  char line[100000];
  assert(fin);
  if(!commandFile) printf("Enter a command: ");
  while(fgets(line, 99999, fin)) {
    SendJSONCommand(line);
    if(!commandFile) printf("Enter a command: ");
  }
  
  if(commandFile)
    fclose(fin);

  networking::cleanup();


  return EXIT_SUCCESS;
}


void SendJSONCommand(std::string queryStr) {
    std::string responseStr;
    
    Json::Rpc::TcpClient tcpClient(std::string("127.0.0.1"), 8086);
    Json::Value query;
    if(!tcpClient.Connect()) {
      std::cerr << "Cannot connect to remote peer!" << std::endl;
      exit(EXIT_FAILURE);
    } 

    std::cout << "Query is: " << queryStr << std::endl;

    if(tcpClient.Send(queryStr) == -1) {
      std::cerr << "Error while sending data!" << std::endl;
      exit(EXIT_FAILURE);
    }

    /* wait the response */
    if(tcpClient.Recv(responseStr) != -1) {
      std::cout << "Received: " << responseStr << std::endl;
    } else {
      std::cerr << "Error while receiving data!" << std::endl;
    }
    
    tcpClient.Close();
}