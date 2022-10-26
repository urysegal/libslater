import sys
from sage.all import *

L_MAX = 8
M_COUNT = ((L_MAX*2) + 1)

with open("gaunt-table.h.in") as f:
    lines = f.readlines()
    with open("gaunt-table.h", "w") as r:
        r.write("#pragma once\n")
        r.write(f"#define L_MAX {L_MAX}\n")
        for line in lines:
            r.write(line)

with open("gaunt-table.cpp.in") as f:
    lines = f.readlines()
    with open("gaunt-table.cpp", "w") as r:
        for line in lines:
            r.write(line)
        done=0

        for l1 in range(0,L_MAX+1):
            r.write("\t{\n")
            for l2 in range(0,L_MAX+1):
                r.write("\t\t{\n")
                for l3 in range(0,L_MAX+1):
                    r.write("\t\t\t{\n")
                    for m1 in range(-L_MAX,L_MAX+1):
                        r.write("\t\t\t\t{\n")
                        for m2 in range(-L_MAX,L_MAX+1):
                            r.write(f"\t\t\t\t\t{{ /* l1={l1} l2={l2} l3={l3} m1={m1} m2={m2} */ ")
                            for m3 in range(-L_MAX,L_MAX+1):
                                x = float(gaunt(l1,l2,l3,m1,m2,m3))
                                r.write(f"{x}, ")
                                done = done + 1
                                if done % 100000 == 0 :
                                    print(f"{done} done\n")
                            r.write("\t\t\t\t\t },\n")
                        r.write("\t\t\t\t},\n")
                    r.write("\t\t\t},\n")
                r.write("\t\t},\n")
            r.write("\t},\n")
        r.write("}; // array gaunt_table \n")
        r.write("} // namespace ")
