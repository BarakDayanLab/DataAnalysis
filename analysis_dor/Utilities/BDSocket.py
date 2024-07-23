import os
import json
import socket
import struct
import threading
from Utilities.BDLogger import BDLogger


# ---------------------------------------------------------------------#
# IMPORTANT NOTE
#
# If you wish to enable another connection, not that it must be with
# the 5050-5059 ports. This is what the FireWall in lab room 112
# is configured with. You can change the rule in Windows Defender Firewall
# ---------------------------------------------------------------------#

class BDSocket:

    SERVER_IP_DROR = "132.77.55.172"  # Dror's machine at room 113
    SERVER_IP_LAB = "132.77.54.212"  # Lab's machine at room 112
    SERVER_IP_LOOPBACK = "127.0.0.1"  # Localhost

    PORT = 5050  # server port number

    def __init__(self, connections_map, writeable, server_id):

        self.logger = BDLogger()

        self.server_id = server_id

        self.writeable = writeable

        self.hostname = socket.gethostname()
        self.my_ip_address = socket.gethostbyname(self.hostname)

        self.connections_map = connections_map
        pass

    def handle_client(self, client_socket, addr, connection_name, writeable):
        try:
            # Get the relevant entry
            connection = self.connections_map[connection_name]

            request = None
            request_unpacked = None

            while True:

                # Receive client messages
                data_size = client_socket.recv(4)
                data_size = struct.unpack("!I", data_size)[0]
                request = client_socket.recv(data_size).decode("utf-8")  # 1024
                request_object = json.loads(request)

                if 'callback' in connection:
                    connection['callback'](request_object)
                    continue
                else:
                    # Assign the value received
                    writeable[connection['value_name']] = request_object

                # Send a message back to the client
                response_message = f"accepted:{self.server_id}"
                response_message = json.dumps(response_message)
                response_message = response_message.encode('utf-8')
                data_size = len(response_message)
                struct_data = struct.pack('!I', data_size) + response_message
                client_socket.send(struct_data)

        except Exception as e:
            self.logger.error(f"Error when handling client: {e}.")
            connection['thread'] = None
        finally:
            client_socket.close()
            self.logger.info(f"Connection to client ({addr[0]}:{addr[1]}) closed")

    def run_server(self):
        self.logger.info(f'Socket listening on host {self.hostname} on {self.my_ip_address}:{BDSocket.PORT}')
        server_thread = threading.Thread(target=self._run_server)
        server_thread.start()

    def _get_connection_entry(self, address, port):
        """ Look for the connection that has this address defined
            Return None if not found (e.g. - not an approved connection)
        """
        # Iterate over all connections
        for attribute, entry in self.connections_map.items():
            if '*' in entry['trusted_IP']:
                return attribute, entry
            if address in entry['trusted_IP']:
                return attribute, entry
        return None, None

    def _run_server(self):
        # Create a socket object
        server = None
        try:
            server = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
            # bind the socket to the host and port
            server.bind((self.my_ip_address, BDSocket.PORT))
            # listen for incoming connections
            server.listen()
            self.logger.info(f'Listening on {self.my_ip_address}:{BDSocket.PORT}')

            while True:
                # accept a client connection
                client_socket, addr = server.accept()

                # Check if this is a connection we are ready to accept
                connection_name, connection_entry = self._get_connection_entry(addr[0], addr[1])
                if connection_entry is None:
                    self.logger.warn(f'Attempt to connect from {addr[0]}:{addr[1]}. Not approved!')
                    continue

                # Accept the connection & start the new thread to handle it
                self.logger.info(f"Accepted connection from {addr[0]}:{addr[1]}")
                connection_entry['thread'] = threading.Thread(target=self.handle_client, args=(client_socket, addr, connection_name, self.writeable))
                connection_entry['thread'].start()

        except Exception as e:
            self.logger.error(e)
        finally:
            server.close()

    @staticmethod
    def run_client():
        """
        Run a test client - connect to the server and send a message
        """
        # Create a socket object
        client = socket.socket(socket.AF_INET, socket.SOCK_STREAM)

        # Establish connection with server
        client.connect((BDSocket.SERVER_IP, BDSocket.PORT))

        try:
            while True:
                # get input message from user and send it to the server
                msg = input("Enter message: ")
                client.send(msg.encode("utf-8")[:1024])

                # receive message from the server
                response = client.recv(1024)
                response = response.decode("utf-8")

                # If server sent us "closed" in the payload, we break out of
                # the loop and close our socket
                if response.lower() == "closed":
                    break

                print(f"Received: {response}")
        except Exception as e:
            print(f"Error: {e}")
        finally:
            # close client socket (connection to the server)
            client.close()
            print("Connection to server closed")

if __name__ == "__main__":

    run_as = "client"  # "server"

    if run_as == "server":
        writeable = {}
        bdsocket = BDSocket(writeable)
        bdsocket.run_server()
    else:
        # Run a test client
        BDSocket.run_client()

    pass
