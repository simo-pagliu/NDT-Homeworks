import numpy as np
import matplotlib.pyplot as plt
import numpy as np

with open("fission_CS.txt", "r") as f:
    cons = f.readlines()

text_cs = []
local_cs = []
for con in cons:
    if con == "\n":
        text_cs.append(local_cs)
        local_cs = []
    else:
        local_cs.append(con.split())

fission_cs = []
for cs in text_cs:
    print(cs[0][0])
    cs_dict = {
        "Z": int(cs[0][0]),
        "A": int(cs[0][1]),
        "E": [],
        "s1": [],
        "s2": [],
    }
    for i in range(1, len(cs)):
        cs_dict["E"].append(float(cs[i][0]))
        cs_dict["s1"].append(float(cs[i][1]))
        cs_dict["s2"].append(float(cs[i][2]))

    fission_cs.append(cs_dict)

def get_fission_cs(Z, A):
    for cs in fission_cs:
        if cs["Z"] == Z and cs["A"] == A:
            return cs
    return None

for ZA in ((92, 235), (92, 238), (94, 239)):
    isotope = get_fission_cs(*ZA)
    plt.loglog(isotope["E"], isotope["s1"], label=f"{ZA[0]}-{ZA[1]}")
plt.legend()
plt.show()
