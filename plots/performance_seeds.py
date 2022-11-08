import matplotlib.pyplot as plt
import pandas as pd
import os
import numpy as np

full = pd.read_csv("perf_results.csv", index_col="disease")

paths = [filename for filename in os.listdir("../data/diseases")]
seed_counts = []
for p in paths:
  f = open("../data/diseases/" + p, "r")
  seed_counts.append(len(f.readlines()))
  f.close()

full["seed_counts"] = seed_counts

full["average"] = (full["dm_algo"] + full["diamond"] + full["best"] + full["random"])/4

print(full["seed_counts"].corr(full["average"]))

full.plot.scatter(x = "seed_counts", y = "average", s = 1)
plt.title("Seed Counts and Average Performance")
plt.xlabel("Seed Counts")
plt.ylabel("Average Performance")

z = np.polyfit(x=full["seed_counts"], y=full["average"], deg=1)
p = np.poly1d(z)
full['trendline'] = p(full["seed_counts"])

full.set_index("seed_counts", inplace=True)
full.trendline.sort_index(ascending=False).plot()


plt.show()


print(list(full["seed_counts"]).corr(list(full["average"])))