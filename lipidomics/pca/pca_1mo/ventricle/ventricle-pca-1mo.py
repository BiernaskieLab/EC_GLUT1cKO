__author__ = "James Colter"
__copyright__ = "Open Source"
__credits__ = ["James Colter"]
__license__ = "GPL"
__version__ = "1"
__maintainer__ = "James Colter"
__email__ = "jdcolter@ucalgary.ca"
__status__ = "Research"

import matplotlib.pyplot as plt
import matplotlib.ticker as mtick
import seaborn as sns
import numpy as np
import pandas as pd

sns.set_theme(style='ticks', font_scale=2)

data = pd.read_csv(r'ventricle-pca-1mo-Scores.csv', delimiter=';')

ko = data.iloc[:281,:]
ctrl = data.iloc[281:,:]

fig = plt.figure()
plt.rcParams['figure.dpi']=300
plt.rc('axes', labelsize=12)
ax = fig.add_subplot(projection='3d')
ax.scatter(ctrl.iloc[:,1], ctrl.iloc[:,2], ctrl.iloc[:,3], color='darkgray', s=8, marker='o', label=None)
ax.scatter(ko.iloc[:,1], ko.iloc[:,2], ko.iloc[:,3], color='darkred', s=8, marker='^', label=None)
ax.scatter([0],[0], color='darkred', s=80, marker='^', label='Glut1KO')
ax.scatter([0],[0], color='darkgray', s=80, marker='o', label='Control')
ax.set_zlabel('Component 3')
ax.set_ylabel('Component 2')
ax.set_xlabel('Component 1')
plt.legend(loc='upper left')
plt.tight_layout()
ax.xaxis.set_pane_color((1.0, 1.0, 1.0, 1.0))
ax.yaxis.set_pane_color((1.0, 1.0, 1.0, 1.0))
ax.zaxis.set_pane_color((1.0, 1.0, 1.0, 1.0))

plt.show()