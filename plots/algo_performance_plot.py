import matplotlib.pyplot as plt
import pandas as pd

df = pd.read_csv("algo_perf.csv")


f = plt.figure(1, figsize=(50,20))
ax = f.add_subplot(1,1,1)
df.plot(kind='bar', ax = ax)
plt.ylabel("Average Prediction Accuracy")
plt.xlabel("Diseases")
labels = ["Alzheimer's Disease", "Multiple Sclerosis", "Rheumatoid Arthiritis", "Blood Coagulation Disorders", "Blood Platelet Disorders"]
ax.set_xticklabels(labels, rotation=0, fontsize=6)
plt.title("Performance of Disease Module finding methods")
plt.show()