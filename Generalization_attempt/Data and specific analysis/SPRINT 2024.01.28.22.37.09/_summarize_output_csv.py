import pandas as pd

df = pd.read_csv('output_total.csv')

df_summary = df[df['num_data_points'] > 0]
df_summary = df_summary.to_csv('output_summary.csv', index=False)
