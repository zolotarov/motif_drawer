import matplotlib.pyplot as plt

# used to set minor ticks on the x-axis
from matplotlib.ticker import AutoMinorLocator

from Bio import SeqIO
from sys import argv
import re
from collections import defaultdict
from collections import OrderedDict
from operator import itemgetter

# load motifs from a file into a list
file_in = open(argv[1], 'r')
motifs = [line.strip() for line in open(argv[2], 'r')]


# create a name for the chart file output
figure_output_name = "%s.png" % (file_in.name.split('.')[0])

sequences = SeqIO.parse(file_in, 'fasta')

seq_info = list()
motifs_dict = defaultdict(list)
motifs_revcomp_dict = defaultdict(list)
complement_dict = {"A": "T", "T": "A", "C": "G", "G": "C"}


def rev_comp(sequence):
    '''
    reverse complement a sequence
    '''
    return "".join([complement_dict[i] for i in sequence][::-1])

count = 0
for seq in sequences:
    seq_info.append({'length': len(str(seq.seq)), 'id': str(seq.id),
                     'number': count})
    count += 1
    for motif in motifs:
        motifs_dict[motif].append([
            -len(str(seq.seq)) + m.start()
            for m in re.finditer(motif, str(seq.seq))])
        motifs_revcomp_dict[motif].append([
            -len(str(seq.seq)) + m.start()
            for m in re.finditer(rev_comp(motif), str(seq.seq))])
number_of_lines = seq_info[-1]['number'] + 1
max_length = max([seq['length'] for seq in seq_info])


def remove_empty_motifs(motifs_dictionary):
    '''
    remove motifs without any matches
    '''
    motifs_dict_copy = motifs_dictionary.copy()
    for motif in motifs_dictionary:
        empty_count = [1 for i in motifs_dictionary[motif] if not i]
        if sum(empty_count) == number_of_lines:
            motifs_dict_copy.pop(motif, None)
    motifs_dict = motifs_dict_copy
    return motifs_dict

motifs_dict = remove_empty_motifs(motifs_dict)
motifs_revcomp_dict = remove_empty_motifs(motifs_revcomp_dict)

# setup the plot figure, here the size of the output can be changed
fig = plt.figure(1, figsize=[15, 9], facecolor='w', edgecolor='k',
                 frameon=False)
ax = fig.add_subplot(111)

# the axes depend on the maximum sequence length (x)
# and the number of sequences (y)
plt.axis([0, -max_length, 0, number_of_lines+2], frameon=False)
ax.xaxis.grid()
ax.spines['top'].set_visible(False)
ax.spines['left'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.xaxis.set_minor_locator(AutoMinorLocator(4))
plt.tick_params(
    axis='y',          # changes apply to the y-axis
    which='both',      # both major and minor ticks are affected
    left='off',      # ticks along the bottom edge are off
    right='off',         # ticks along the top edge are off
    labelbottom='off',  # labels along the bottom edge are off
    pad=10)
plt.tick_params(
    axis='x',
    which='both',
    bottom='on',
    top='off',
    labelbottom='on')

# start ticks at 1, not 0 and label them according to the deflines (id)
# of the sequences
plt.yticks(range(1, number_of_lines+1), [seq['id'] for seq in seq_info])

# invert the x axis so that it goes from a negative number to 0 on the right
plt.gca().invert_xaxis()


def draw_line(number, length):
    # (y location, x start, x end, ...)
    plt.hlines(number+1, 0, -length, "black", linewidth=3, zorder=0)

for line in seq_info:
    draw_line(line['number'], line['length'])


def draw_motifs(motif, line_number, locations, color, marker):
    plt.scatter(locations, [line_number+1]*len(locations), s=150, c=color,
                marker=marker, alpha=1, edgecolors='k', linewidth=0.3,
                label=motif)


colors = {0: '#99A0A9', 1: '#FF6B6B', 2: '#C7F464', 3: '#a79d88',
          4: '#BD1550', 5: '#E97F02', 6: '#F8CA00', 7: '#8A9B0F',
          8: '#AEE239', 9: '#F307BE', 10: '#1784D6', 11: '#FCB84E',
          12: '#BE9C6F'}

motif_numbers = {v: k for (k, v) in enumerate(motifs_dict.keys())}

for motif in motifs_revcomp_dict:
    for m in enumerate(motifs_revcomp_dict[motif]):
        draw_motifs(motif, m[0], m[1], colors[motif_numbers[motif]], "<")
for motif in motifs_dict:
    for m in enumerate(motifs_dict[motif]):
        draw_motifs(motif, m[0], m[1], colors[motif_numbers[motif]], ">")

plt.xlabel('Distance from start codon')
handles, labels = plt.gca().get_legend_handles_labels()
by_label = OrderedDict(zip(labels, handles))
plt.legend(by_label.values(), by_label.keys(), loc='upper center',
           ncol=5, scatterpoints=1, scatteryoffsets=[0.65]*len(motifs))
plt.savefig(figure_output_name, format="png")
print "Saved the figure in %s" % figure_output_name
# plt.show()
file_in.close()
