#!/usr/bin/env python3

import math
import sys

class point(tuple):
    def __add__(self, o):
        return point([a+b for a,b in zip(self, o)])
    def __sub__(self, o):
        return point([a-b for a,b in zip(self, o)])
    def __mul__(self, s):
        return point([a*s for a in self])
    def __str__(self):
        return "("+",".join([f"{a:.6f}" for a in self])+")"

def angle(d):
    return point((math.cos(d * math.pi / 180), math.sin(d * math.pi / 180)))

x = point((-27.506334, 153.17907))
R = 0.01

def printrow(p,name):
    print(f'"{name}",{p[0]:.6f},{p[1]:.6f},"",0')

def subtri(tri,name,depth):
    (a,b,c)=tri
    m = (a+b+c)*(1/3)
    printrow(m,name)
    if depth > 0:
        subtri((m,a,b),name+"a",depth-1)
        subtri((m,b,c),name+"b",depth-1)
        subtri((m,c,a),name+"c",depth-1)

if __name__ == '__main__':
    d = int(sys.argv[1])
    tri=(x + angle(0)*R,
         x + angle(360/3)*R,
         x + angle(2*360/3)*R)
    print("Name, Latitude, Longitude")
    printrow(tri[0],"A")
    printrow(tri[1],"B")
    printrow(tri[2],"C")
    if d > 1:
        subtri(tri, "o", d - 2)
