import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
sns.set()

# Load the example flights dataset and convert to long-form
normal_data = np.random.randn(10, 12)
ax = sns.heatmap(normal_data, center=0,annot=True)
plt.show()