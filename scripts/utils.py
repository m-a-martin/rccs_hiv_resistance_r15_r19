import numpy as np


def plot_style(grey='#333333'):
    import matplotlib as mpl
    mpl.rcParams['font.family'] = 'sans-serif'
    mpl.rcParams['font.sans-serif'] = 'helvetica'
    mpl.rcParams['font.weight'] = 'light'
    mpl.rcParams['text.color'] = grey
    mpl.rcParams['axes.labelcolor'] = grey
    mpl.rcParams['xtick.color'] = '#707070'
    mpl.rcParams['ytick.color'] = '#707070'
    # Font sizes
    mpl.rcParams['figure.titlesize'] = 16
    mpl.rcParams['axes.titlesize'] = 16
    mpl.rcParams['axes.labelsize'] = 16
    mpl.rcParams['xtick.labelsize'] = 16
    mpl.rcParams['ytick.labelsize'] = 16
    # Border colors
    mpl.rcParams['axes.edgecolor'] = grey
    mpl.rcParams['grid.color'] = '#eaeaea'
    #mpl.rcParams['grid.zorder'] = 0
    # Legend
    mpl.rcParams['legend.fontsize'] = 16
    mpl.rcParams['legend.frameon'] = False
    
def read_fasta(fp):
    name, seq = None, []
    for line in fp:
        line = line.rstrip()
        if line.startswith(">"):
            if name: yield (name, ''.join(seq))
            name, seq = line[1:], []
        else:
            seq.append(line)
    if name: yield (name, ''.join(seq))


def import_aln(fh):
    s_names = []
    all_s = ''
    #fh = open(fasta_path, 'rt')
    with fh as fasta:
        for h,s in read_fasta(fasta):
            s_names.append(h)
            all_s += s
    fh.close()
    n_seqs = len(s_names)
    size = int(len(all_s)/n_seqs)
    s_arr = np.array(list(all_s.lower())).reshape((n_seqs, size))
    return(np.array(s_names), s_arr)


def import_seqs(fasta_path):
    s_names = []
    s_list = []
    fh = open(fasta_path, 'rt')
    with fh as fasta:
        for h,s in read_fasta(fasta):
            s_names.append(h)
            s_list.append(s)
    fh.close()    
    return(np.array(s_names), s_list)

