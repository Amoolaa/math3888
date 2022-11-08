import matplotlib.pyplot as plt
import pandas as pd

full = pd.read_csv("perf_results.csv", index_col="disease")

print(full)

selected = ["blood coagulation disorders", "blood platelet disorders", "lipid metabolism disorders", "nutritional and metabolic diseases", "crohn disease"]

df = full.loc[selected]

print(df)


f = plt.figure(1, figsize=(50,20))
ax = f.add_subplot(1,1,1)
df.plot(kind='bar', ax = ax)
plt.ylabel("Average Prediction Accuracy")
plt.xlabel("Diseases")
ax.set_xticklabels(selected, rotation=0, fontsize=6)
plt.legend(["DIAMOnD", "PCSTP-DM", "Best Neighbours", "Random Neighbours"])
plt.title("Performance of Disease Module finding methods")
plt.show()