import numpy as np 
import socket, struct, time

try:
    import matplotlib.pylab as plt
except:
    print("Couldn't import matplotlib")

server_socket = socket.socket(socket.AF_INET, socket.SOCK_DGRAM)
server_socket.bind(("", 31664))

def listen(fn=None):
    if fn is not None:
        f = open(fn, 'w+')

    while True:
        recvdata, address = server_socket.recvfrom(2400)
        print(address, recvdata[:88])
        if fn is not None:
            f.write(recvdata)
