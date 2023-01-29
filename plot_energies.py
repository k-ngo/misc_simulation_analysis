import matplotlib
import matplotlib.pylab as plt
from matplotlib.ticker import AutoMinorLocator
import seaborn as sns
import glob
import pandas as pd
import os
sns.set(style='ticks', context='talk', palette='Set1')
colors = [sns.color_palette('Set1')[0], sns.color_palette('Set1')[1], sns.color_palette('Set1')[2], sns.color_palette('Set1')[3], sns.color_palette('Set1')[4], sns.color_palette('Set1')[5]]
colors_energy = [sns.color_palette('Set1')[0], 'maroon', ]

fig, axes = plt.subplots(4, 2, figsize=(12, 14))

for folder, label, ax, num in zip(['deactVSD2', 'actVSD2', 'deactVSD4', 'actVSD4'], ['Deact. VSD II', 'Act. VSD II', 'Deact. VSD IV', 'Act. VSD IV'], axes.flatten(), range(4)):
    time = 1000
    print(num)

    solvation_free_energy = pd.read_csv(os.path.join(folder, 'withmem', 'en.final.txt'), header=None)
    interaction_energy_in_vacuum = pd.read_csv(os.path.join(folder, 'withoutmem', 'entro.final.txt'), header=None)
    RMS_complex = pd.read_csv(os.path.join(folder, 'rmsd_COMPLEX.txt'), header=0, sep='\s+').fillna(0)
    RMS_TOX = pd.read_csv(os.path.join(folder, 'rmsd_TOX.txt'), header=0, sep='\s+').fillna(0)
    RMS_VSD = pd.read_csv(os.path.join(folder, 'rmsd_VSD.txt'), header=0, sep='\s+').fillna(0)

    # Convert frame to time
    solvation_free_energy.index = [round(i / solvation_free_energy.shape[0] * time, 2) for i in solvation_free_energy.index.values]
    interaction_energy_in_vacuum.index = [round(i / interaction_energy_in_vacuum.shape[0] * time, 2) for i in interaction_energy_in_vacuum.index.values]
    RMS_complex.index = [round(i / RMS_complex.shape[0] * time, 2) for i in RMS_complex.index.values]
    RMS_TOX.index = [round(i / RMS_TOX.shape[0] * time, 2) for i in RMS_TOX.index.values]
    RMS_VSD.index = [round(i / RMS_VSD.shape[0] * time, 2) for i in RMS_VSD.index.values]

    # Plot raw data
    # sns.lineplot(data=RMS_complex['mol0'], ax=ax, color=colors[0], alpha=0.2)
    sns.lineplot(data=RMS_TOX['mol0'], ax=ax, color='mediumseagreen', alpha=1, label='PTx2')
    sns.lineplot(data=RMS_VSD['mol0'], ax=ax, color='dimgray', alpha=1, label=label)

    if num == 0 or num == 1:
        color = colors[0] if num == 0 else colors[1]
        sns.lineplot(data=solvation_free_energy[0], ax=axes[2, 0], color=color, alpha=0.2, legend=False)
        sns.lineplot(data=interaction_energy_in_vacuum[0], ax=axes[3, 0], color=color, alpha=0.2, legend=False)
    else:
        color = colors[0] if num == 2 else colors[1]
        sns.lineplot(data=solvation_free_energy[0], ax=axes[2, 1], color=color, alpha=0.2, legend=False)
        sns.lineplot(data=interaction_energy_in_vacuum[0], ax=axes[3, 1], color=color, alpha=0.2, legend=False)

    # Rolling average calculations
    solvation_free_energy[0] = solvation_free_energy[0].rolling(5).mean()
    interaction_energy_in_vacuum[0] = interaction_energy_in_vacuum[0].rolling(5).mean()
    RMS_complex['mol0'] = RMS_complex['mol0'].rolling(5).mean()
    RMS_TOX['mol0'] = RMS_TOX['mol0'].rolling(5).mean()
    RMS_VSD['mol0'] = RMS_VSD['mol0'].rolling(5).mean()

    # Plot rolling average
    # sns.lineplot(data=RMS_complex['mol0'], ax=ax, color=colors[0], label='Complex')
    # sns.lineplot(data=RMS_TOX['mol0'], ax=ax, color=colors[2], label='PTx2')
    # sns.lineplot(data=RMS_VSD['mol0'], ax=ax, color='slategray', label=label)
    if num == 0 or num == 1:
        color = colors[0] if num == 0 else colors[1]
        sns.lineplot(data=solvation_free_energy[0], ax=axes[2, 0], color=color, label=label)
        sns.lineplot(data=interaction_energy_in_vacuum[0], ax=axes[3, 0], color=color, label=label)
    else:
        color = colors[0] if num == 2 else colors[1]
        sns.lineplot(data=solvation_free_energy[0], ax=axes[2, 1], color=color, label=label)
        sns.lineplot(data=interaction_energy_in_vacuum[0], ax=axes[3, 1], color=color, label=label)

    # axes[num, 0].set_title('Distribution of Interactions', fontsize=18, y=1.04)
    # axes[num, 0].axis('off')

    # axes[num, 1].set_title('RMSD from Starting Structure', fontsize=18, y=1.04)
    ax.set_xlabel('Time (ns)', fontsize=20)
    ax.set_ylabel('RMSD (Å)', fontsize=20)
    ax.legend(frameon=False, fontsize=13)  # , bbox_to_anchor=(0, 1.02, 1, 0.2), loc='upper left', ncol=2)
    ax.set_ylim(-0.6, 5.3)

# Labels
# axes[4, 0].set_title('MM/PBSA Enthalpy', fontsize=18, y=1.04)
axes[2, 0].set_xlabel('Time (ns)', fontsize=20)
axes[2, 0].set_ylabel('Enthalpy (kcal/mol)', fontsize=20)
axes[2, 0].legend().remove()
axes[2, 0].legend(frameon=False, fontsize=13, bbox_to_anchor=(0, 1.02, 1, 0.2), loc='upper left', ncol=2)
axes[2, 0].set_ylim(-79, 45)

axes[3, 0].set_xlabel('Time (ns)', fontsize=20)
axes[3, 0].set_ylabel('Interaction Energy\n(kcal/mol)', fontsize=20)
axes[3, 0].legend().remove()
axes[3, 0].legend(frameon=False, fontsize=13, bbox_to_anchor=(0, 1.02, 1, 0.2), loc='upper left', ncol=2)
axes[3, 0].set_ylim(-220, -45)

# axes[4, 1].set_title('Interaction Energy in Vacuum', fontsize=18, y=1.04)
axes[2, 1].set_xlabel('Time (ns)', fontsize=20)
axes[2, 1].set_ylabel('Enthalpy (kcal/mol)', fontsize=20)
# axes[2, 1].legend(frameon=False, fontsize=13, bbox_to_anchor=(1.01, 1), loc='upper left')
axes[2, 1].legend(frameon=False, fontsize=13, bbox_to_anchor=(0, 1.02, 1, 0.2), loc='upper left', ncol=2)
axes[2, 1].set_ylim(-79, 45)

axes[3, 1].set_xlabel('Time (ns)', fontsize=20)
axes[3, 1].set_ylabel('Interaction Energy\n(kcal/mol)', fontsize=20)
# axes[3, 1].legend(frameon=False, fontsize=13, bbox_to_anchor=(1.01, 1), loc='upper left')
axes[3, 1].legend(frameon=False, fontsize=13, bbox_to_anchor=(0, 1.02, 1, 0.2), loc='upper left', ncol=2)
axes[3, 1].set_ylim(-220, -45)

for ax in axes.flatten():
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.yaxis.set_minor_locator(AutoMinorLocator())
    ax.axvspan(0, 90, facecolor='mediumseagreen', alpha=.2)

sns.despine()
plt.tight_layout()
plt.savefig('figure.png', dpi=300)
