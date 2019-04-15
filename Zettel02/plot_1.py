import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd

plt.rc('text', usetex=True)
plt.rc('font', family='serif')

df = pd.read_csv('Aufgabe1/bild_10.txt', decimal='.', delimiter=';')
print(df.head())

# beispeil für heatmap aus dem Internet
Index = ['aaa', 'bbb', 'ccc', 'ddd', 'eee']
Cols = ['A', 'B', 'C', 'D']
df = pd.DataFrame(abs(np.random.randn(5, 4)), index=Index, columns=Cols)
print(df.head())
#array = np.genfromtxt("Aufgabe1/bild_10.txt", unpack=True)
#print(array)

#ax = sns.heatmap(df, annot=True)
#plt.show()
