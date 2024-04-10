import numpy as np
import pandas as pd

df_experiment_summary = pd.read_csv('output_experiment_summary.csv')
df_normalization_summary = pd.read_csv('output_normalization_summary.csv')
N_cycles_experiment = 250
N_cycles_normalization = 50

SNR = (df_experiment_summary['experiment_counts'][(x.endswith('transit') for x in df_experiment_summary['x_vals'])] / N_cycles_experiment) / (
            df_normalization_summary['experiment_counts'][(x.endswith('transit') for x in df_normalization_summary['x_vals'])] / N_cycles_normalization)

# print(df_experiment_summary)
transmissions = (df_experiment_summary['experiment_counts'][(x.endswith('transmission') for x in df_experiment_summary['x_vals'])])
reflections = (df_experiment_summary['experiment_counts'][(x.endswith('reflection') for x in df_experiment_summary['x_vals'])])

fidelities_T = np.array(transmissions) / (np.array(transmissions) + np.array(reflections))

print(fidelities_T)
print(np.array(SNR))

# print(SNR)
