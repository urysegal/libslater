L_MAX = 10


with open("gaunt-table.h.in") as f:
    lines = f.readlines()
    with open("gaunt-table.h", "w") as r:
        r.write("#pragma once\n")
        r.write(f"#define L_MAX {L_MAX}\n")
        for line in lines:
            r.write(line)

