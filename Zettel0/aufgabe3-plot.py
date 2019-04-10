import matplotlib.pyplot as plt
import pandas as pd

df = pd.read_table('euler_vergleich.txt', skiprows=1, decimal='.', delimiter=';')
print(df.head())
