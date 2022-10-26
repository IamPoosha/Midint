from collections import defaultdict
import numpy as np
import math
import matplotlib
import matplotlib.pyplot as plt
import os
import re
import subprocess
import time

# before starting, make sure "main.cpp" is compiled to "wide"
time_each_benchmark = 1
start_time = time.time()
output = subprocess.run(['./wide', str(time_each_benchmark)], stdout=subprocess.PIPE)
end_time = time.time()
lines = output.stdout.decode("utf-8").split("\n")
# TODO: update regex!
pattern = r"ops per second:  (\w+)~(\w+)~(\d+) : 2\^ ([^\s]+)    \|\|    time: ([^\s]+)"
regex = re.compile(pattern)
points = defaultdict(lambda: defaultdict(list))

for line in lines:
    match = regex.findall(line)
    if match:
        test_name, int_type, operand_size, ops_ps, test_time = match[0]
        operand_size = int(operand_size)
        ops_ps = float(ops_ps)
        test_time = float(test_time)
        points[test_name][int_type].append((operand_size, ops_ps, test_time))

for (test_name, int_types) in points.items():
    plt.cla()
    plt.figure(figsize=(14, 8), dpi=600)
    # having x in x scale, the poor way.
    xticks = [64,128,192,256,384,512,768,1024]
    logged_xticks = [math.log(tick) for tick in xticks]
    for (int_type, pts) in int_types.items():
        pts.sort()
        if not all(0.6 < test_time / time_each_benchmark < 1.6
                   for operand_size, ops_ps, test_time in pts
        ):
            # 31.1 is log2(get_num_cycles_per_second(std::exp2(30)))
            rates = [round(ops_ps-32.2, 2)
                     for operand_size, ops_ps, test_time in pts]
            print(f"adjust rates for {test_name}~{int_type}",
                  str(rates).replace("[", "{").replace("]", "}"))
        xs = [p[0] for p in pts]
        assert xs == xticks
        ys = [p[1] for p in pts]
        plt.plot(logged_xticks, ys,label=f"{test_name}~{int_type}",linewidth=3)
    plt.legend()
    plt.xticks(logged_xticks, xticks)
    start,end = plt.gca().get_ylim()
    plt.gca().yaxis.set_ticks(np.arange(round(start), round(end), 1))
    plt.gca().set_xlabel('operand size (bits)')
    plt.gca().set_ylabel('lg2 ops/s')
    plt.grid()
    file_name = f"./figures/{test_name}.png"
    plt.savefig(file_name)

print("time factor", (end_time-start_time)/time_each_benchmark)
