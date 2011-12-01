function str = tcp_request(host, port, req)
    import java.io.*
    import java.net.*;

    sock = [];
    str = '';

    try
        sock = Socket(host, port);
        in = BufferedReader(InputStreamReader(sock.getInputStream));
        out = PrintWriter(sock.getOutputStream,true);
        out.println(req);
        while true 
            line=in.readLine();
            if isempty(line), break; end
            str = strcat(str, char(line));
        end
        out.close();
        in.close();
        sock.close();
    catch e
        e
        e.message
        if ~isempty(sock)
            sock.close();
        end
    end
end