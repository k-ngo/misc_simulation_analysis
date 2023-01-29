import matplotlib
import matplotlib.pylab as plt
from matplotlib.ticker import AutoMinorLocator
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import numpy as np
import seaborn as sns
import glob
import pandas as pd
from matplotlib.axes import Axes
from matplotlib.figure import Figure
import os
sns.set(style='ticks', context='talk', palette='Set1')

skip = 45
# cutoff = 8
cutoff_PTx2 = 10
cutoff_VSD = 5
ligand_folder = 'ligand_all'
receptor_folder = 'receptor'
output = 'energy'
time_total = 1000
col_names = ['frame', 'elec', 'vdw']
# col_names = ['frame', 'bond', 'angle', 'dih', 'vdw14', 'elec14', 'vdw', 'elec', 'total']


def errplot(x, y, yerr, **kwargs):
    ax = plt.gca()
    data = kwargs.pop('data')
    data.plot(x=x, y=y, yerr=yerr, rot=0, kind='bar', ax=ax, alpha=0, **kwargs)


fig1, axes1 = plt.subplots(1, 1, figsize=(10, 5), num='mean')

# Read ligand data
ligand = pd.DataFrame()
for sim, label, num in zip(['deactVSD2', 'actVSD2', 'deactVSD4', 'actVSD4'],
                           ['Deact.\nVSD II', 'Act.\nVSD II', 'Deact.\nVSD IV', 'Act.\nVSD IV'],
                           [0, 0.1, 0.2, 0.3]):
    # Loop through each residue energy file
    df = pd.DataFrame()
    df['Residue'] = [os.path.splitext(os.path.basename(file))[0] for file in glob.glob(os.path.join(sim, ligand_folder, '*.dat'))]
    for file in glob.glob(os.path.join(sim, ligand_folder, '*.dat')):
        resid = os.path.splitext(os.path.basename(file))[0]
        data = pd.read_csv(file, header=None, sep='\s+', skiprows=1 + skip, names=col_names)
        df.loc[df['Residue'] == resid, 'ID'] = int(resid[1:]) + num
        df.loc[df['Residue'] == resid, 'System'] = label
        df.loc[df['Residue'] == resid, 'Interaction Energy\nw/ VSD + Lipids\n(kcal/mol)'] = data['elec'].mean() + data['vdw'].mean()
        df.loc[df['Residue'] == resid, 'Std Dev'] = (data['elec'] + data['vdw']).std()
    ligand = pd.concat([ligand, df])

# Sort by residue numbers
ligand.sort_values(by='ID', inplace=True)
# Remove residues with interaction energy less than cutoff
few_interactions = ligand.copy()
few_interactions['Interaction Energy\nw/ VSD + Lipids\n(kcal/mol)'] = few_interactions['Interaction Energy\nw/ VSD + Lipids\n(kcal/mol)'].abs()
few_interactions = few_interactions.groupby('Residue').sum()['Interaction Energy\nw/ VSD + Lipids\n(kcal/mol)'].abs() < cutoff_PTx2
ligand = ligand[~ligand['Residue'].isin(few_interactions[few_interactions == True].index)]
print(ligand.to_string())
cat = sns.catplot(x='System', y='Interaction Energy\nw/ VSD + Lipids\n(kcal/mol)', col='Residue', data=ligand,
                  kind='bar', col_wrap=5, height=3, aspect=1.2, sharex=False)
cat.map_dataframe(errplot, 'System', 'Interaction Energy\nw/ VSD + Lipids\n(kcal/mol)', 'Std Dev')
for ax in cat.axes.flat:
    ax.grid(True, axis='y')
    ax.set(xlabel='')
    ax.set_ylim(-44, 5)
    ax.set_ylabel('Interaction Energy\nw/ VSD + Lipids\n(kcal/mol)', rotation='horizontal', va='center', ha='right', fontsize=16)
    ax.set_xticklabels(ax.get_xticklabels(), fontsize=13)  # rotation=15, ha='right', rotation_mode='anchor'
    ax.set_title(r'$\bf{PTx2}$ ' + ax.get_title().split()[-1], fontsize=16)

plt.tight_layout()
plt.savefig('PTx2.png', dpi=300)

# Read receptor data (VSD 2)
receptor = pd.DataFrame()
for sim, label, num in zip(['deactVSD2', 'actVSD2'],
                           ['Deact.\nVSD II', 'Act.\nVSD II'],
                           [0, 0.1]):
    # Loop through each residue energy file
    df = pd.DataFrame()
    df['Residue'] = [os.path.splitext(os.path.basename(file))[0] for file in glob.glob(os.path.join(sim, receptor_folder, '*.dat'))]
    for file in glob.glob(os.path.join(sim, receptor_folder, '*.dat')):
        resid = os.path.splitext(os.path.basename(file))[0]
        data = pd.read_csv(file, header=None, sep='\s+', skiprows=1 + skip, names=col_names)
        df.loc[df['Residue'] == resid, 'ID'] = int(resid[1:]) + num
        df.loc[df['Residue'] == resid, 'System'] = label
        df.loc[df['Residue'] == resid, 'Interaction Energy\nw/ PTx2 Only\n(kcal/mol)'] = data['elec'].mean() + data['vdw'].mean()
        df.loc[df['Residue'] == resid, 'Std Dev'] = (data['elec'] + data['vdw']).std()
    receptor = pd.concat([receptor, df])

# Sort by residue numbers
receptor.sort_values(by='ID', inplace=True)
# Remove residues with interaction energy less than cutoff
few_interactions = receptor.copy()
few_interactions['Interaction Energy\nw/ PTx2 Only\n(kcal/mol)'] = few_interactions['Interaction Energy\nw/ PTx2 Only\n(kcal/mol)'].abs()
few_interactions = few_interactions.groupby('Residue').sum()['Interaction Energy\nw/ PTx2 Only\n(kcal/mol)'].abs() < cutoff_VSD
receptor = receptor[~receptor['Residue'].isin(few_interactions[few_interactions == True].index)]
print(receptor.to_string())
cat = sns.catplot(x='System', y='Interaction Energy\nw/ PTx2 Only\n(kcal/mol)', col='Residue', data=receptor,
                  kind='bar', col_wrap=5, height=3, aspect=1.2, sharex=False)
cat.map_dataframe(errplot, 'System', 'Interaction Energy\nw/ PTx2 Only\n(kcal/mol)', 'Std Dev')
for ax in cat.axes.flat:
    ax.grid(True, axis='y')
    ax.set(xlabel='')
    ax.set_ylim(-44, 5)
    ax.set_ylabel('Interaction Energy\nw/ PTx2 Only\n(kcal/mol)', rotation='horizontal', va='center', ha='right', fontsize=16)
    ax.set_xticklabels(ax.get_xticklabels(), fontsize=14)
    ax.set_title(r'$\bf{VSD\ II}$ ' + ax.get_title().split()[-1], fontsize=16)

plt.tight_layout()
plt.savefig('VSD2.png', dpi=300)

# Read receptor data (VSD 4)
receptor = pd.DataFrame()
for sim, label, num in zip(['deactVSD4', 'actVSD4'],
                           ['Deact.\nVSD IV', 'Act.\nVSD IV'],
                           [0, 0.1]):
    # Loop through each residue energy file
    df = pd.DataFrame()
    df['Residue'] = [os.path.splitext(os.path.basename(file))[0] for file in glob.glob(os.path.join(sim, receptor_folder, '*.dat'))]
    for file in glob.glob(os.path.join(sim, receptor_folder, '*.dat')):
        resid = os.path.splitext(os.path.basename(file))[0]
        data = pd.read_csv(file, header=None, sep='\s+', skiprows=1 + skip, names=col_names)
        df.loc[df['Residue'] == resid, 'ID'] = int(resid[1:]) + num
        df.loc[df['Residue'] == resid, 'System'] = label
        df.loc[df['Residue'] == resid, 'Interaction Energy\nw/ PTx2 Only\n(kcal/mol)'] = data['elec'].mean() + data['vdw'].mean()
        df.loc[df['Residue'] == resid, 'Std Dev'] = (data['elec'] + data['vdw']).std()
    receptor = pd.concat([receptor, df])

# Sort by residue numbers
receptor.sort_values(by='ID', inplace=True)
# Remove residues with interaction energy less than cutoff
few_interactions = receptor.copy()
few_interactions['Interaction Energy\nw/ PTx2 Only\n(kcal/mol)'] = few_interactions['Interaction Energy\nw/ PTx2 Only\n(kcal/mol)'].abs()
few_interactions = few_interactions.groupby('Residue').sum()['Interaction Energy\nw/ PTx2 Only\n(kcal/mol)'].abs() < cutoff_VSD
receptor = receptor[~receptor['Residue'].isin(few_interactions[few_interactions == True].index)]
print(receptor.to_string())
cat = sns.catplot(x='System', y='Interaction Energy\nw/ PTx2 Only\n(kcal/mol)', col='Residue', data=receptor,
                  kind='bar', col_wrap=5, height=3, aspect=1.2, palette=sns.color_palette([sns.color_palette('Set1')[2], sns.color_palette('Set1')[3]]), sharex=False)
cat.map_dataframe(errplot, 'System', 'Interaction Energy\nw/ PTx2 Only\n(kcal/mol)', 'Std Dev')
for ax in cat.axes.flat:
    ax.grid(True, axis='y')
    ax.set(xlabel='')
    ax.set_ylim(-44, 5)
    ax.set_ylabel('Interaction Energy\nw/ PTx2 Only\n(kcal/mol)', rotation='horizontal', va='center', ha='right', fontsize=16)
    ax.set_xticklabels(ax.get_xticklabels(), fontsize=14)
    ax.set_title(r'$\bf{VSD\ IV}$ ' + ax.get_title().split()[-1], fontsize=16)

plt.tight_layout()
plt.savefig('VSD4.png', dpi=300)

exit()

    # ## Receptor
    # receptor = pd.DataFrame()
    # # Loop through each residue energy file
    # for file in glob.glob(os.path.join(sim, receptor_folder, '*.dat')):
    #     resid = os.path.splitext(os.path.basename(file))[0]
    #     data = pd.read_csv(file, header=None, sep='\s+', skiprows=1 + skip, names=col_names)
    #     for index, row in data.iterrows():
    #         receptor.loc[resid, index] = round(data['elec'][index] + data['vdw'][index], 2)
    #     receptor.loc[resid, 'Electrostatic energy'] = round(data['elec'].mean(), 2)
    #     receptor.loc[resid, 'van der Waals energy'] = round(data['vdw'].mean(), 2)
    #     receptor.loc[resid, 'Interaction Energy'] = round(data['elec'].mean() + data['vdw'].mean(), 2)
    # # Remove residues with interaction energy less than cutoff
    # receptor.drop(receptor[(-cutoff < receptor['Interaction Energy']) & (receptor['Interaction Energy'] < cutoff)].index, inplace=True)
    # # Sort by residue number
    # receptor = receptor.reindex(sorted(receptor.index, key=lambda x: int(x[1:])), axis=0)
    # print(receptor)
    # # Exclude certain columns for plotting
    # receptor_barplot = receptor.loc[:, receptor.columns.isin(['Electrostatic energy', 'van der Waals energy'])]
    # receptor_timeseries = receptor.loc[:, ~receptor.columns.isin(['Electrostatic energy', 'van der Waals energy', 'Interaction Energy'])]
    # # Plot mean energy
    # receptor_barplot.plot(ax=axes1.flatten()[ax2], kind='bar', stacked=True)
    # axes1.flatten()[ax2].set_xlabel(label + ' Residues', fontsize=20)
    # # Plot energy as time series
    # my_map = sns.heatmap(receptor_timeseries,
    #                      xticklabels=False, yticklabels=1, ax=axes3.flatten()[ax1],
    #                      cmap=sns.color_palette('magma', as_cmap=True),
    #                      cbar=True)

# Mean energy plot configs
for ax in axes1.flatten():
    ax.set_ylabel('Energy (kcal/mol)')
    ax.set_xticklabels(ax.get_xticklabels(), rotation=45, fontsize=14)
    ax.yaxis.set_minor_locator(AutoMinorLocator())
    ax.tick_params(which='both', width=2)
    ax.tick_params(which='major', length=7)
    ax.tick_params(which='minor', length=4)
    ax.grid(True, axis='y')

for ax in axes1.flatten()[:-2]:
    ax.legend(frameon=False, fontsize=13, bbox_to_anchor=(0, 1.02, 1, 0.2), loc='upper left')

for ax in axes1.flatten()[2:]:
    ax.get_legend().remove()

sns.despine()
plt.tight_layout()
plt.savefig(output + '.png', dpi=300)
