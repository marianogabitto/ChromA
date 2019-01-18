import classes.data_handle as dh
import numpy as np
import pickle
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages


# #####################################################################
# #####################################################################
# #####################################################################
# COMPUTE RAW SINGLE CELL FILE CORRELATION
def correl(a, b):
    a1 = (a - np.mean(a)) / (np.std(a) * len(a))
    b1 = (b - np.mean(b)) / (np.std(b))

    return np.correlate(a1[:, 0], b1[:, 0])


human_correl = []
human_chromlens = dh.human_lens()
human_chromlens.__delitem__('chrM')
human_chromlens.__delitem__('chrX')
human_chromlens.__delitem__('chrY')
for c_ in human_chromlens:
    print(c_)
    f10000hg = dh.chr_reads(['/Users/mgabitto/Dropbox (Simons Foundation)/2_PeakCalling/ATACseq_Data/10x_scTest/GM12878_mix_mm10/hg19/atac_v1_hg_10k_fragments.tsv'], c_, 1, chromlens[c_])
    f5000hg = dh.chr_reads(['/Users/mgabitto/Dropbox (Simons Foundation)/2_PeakCalling/ATACseq_Data/10x_scTest/GM12878_mix_mm10/hg19-5k/atac_v1_hg_5k_fragments.tsv'], c_, 1, chromlens[c_])
    f1000hg = dh.chr_reads(['/Users/mgabitto/Dropbox (Simons Foundation)/2_PeakCalling/ATACseq_Data/10x_scTest/GM12878_mix_mm10/hg19-1k/atac_v1_hg_1k_fragments.tsv'], c_, 1, chromlens[c_])
    f500hg = dh.chr_reads(['/Users/mgabitto/Dropbox (Simons Foundation)/2_PeakCalling/ATACseq_Data/10x_scTest/GM12878_mix_mm10/hg19-500/atac_v1_hg_500_fragments.tsv'], c_, 1, chromlens[c_])
    human_correl.append([correl(f10000hg, f10000hg), correl(f10000hg, f5000hg), correl(f10000hg, f1000hg), correl(f10000hg, f500hg)])
pickle.dump(human_correl, open('human_raw_correl.pkl', 'wb'))

mouse_correl = {}
mouse_chromlens = dh.mouse_lens()
mouse_chromlens.__delitem__('chrM')
mouse_chromlens.__delitem__('chrX')
mouse_chromlens.__delitem__('chrY')
for c_ in mouse_chromlens:
    print('c_')
    f10000mm = dh.chr_reads(['/Users/mgabitto/Dropbox (Simons Foundation)/2_PeakCalling/ATACseq_Data/10x_scTest/GM12878_mix_mm10/mm10/atac_v1_mm_10k_fragments.tsv'], c_, 1, chromlens[c_])
    f5000mm = dh.chr_reads(['/Users/mgabitto/Dropbox (Simons Foundation)/2_PeakCalling/ATACseq_Data/10x_scTest/GM12878_mix_mm10/mm10-5k/atac_v1_mm_5k_fragments.tsv'], c_, 1, chromlens[c_])
    f1000mm = dh.chr_reads(['/Users/mgabitto/Dropbox (Simons Foundation)/2_PeakCalling/ATACseq_Data/10x_scTest/GM12878_mix_mm10/mm10-1k/atac_v1_mm_1k_fragments.tsv'], c_, 1, chromlens[c_])
    f500mm = dh.chr_reads(['/Users/mgabitto/Dropbox (Simons Foundation)/2_PeakCalling/ATACseq_Data/10x_scTest/GM12878_mix_mm10/mm10-500/atac_v1_mm_500_fragments.tsv'], c_, 1, chromlens[c_])
    mouse_correl[c_] = [1., correl(f10000mm, f5000mm), correl(f10000mm, f1000mm), correl(f10000mm, f500mm)]
pickle.dump(mouse_correl, open('mouse_raw_correl.pkl', 'wb'))


# ###############################################################
# ###############################################################
# ###############################################################
# COMPUTE 10X BED CORRELATION
def bed_file_correlation(f1, f2, f3, f4, chromlens, filt=0):
    f1bed = dh.read_bed(f1)
    f2bed = dh.read_bed(f2)
    f3bed = dh.read_bed(f3)
    f4bed = dh.read_bed(f4)

    cc_correl = {}
    for c_ in chromlens:
        f1cc = np.zeros((chromlens[c_], 1))
        f2cc = np.zeros((chromlens[c_], 1))
        f3cc = np.zeros((chromlens[c_], 1))
        f4cc = np.zeros((chromlens[c_], 1))

        last = 0
        for b_ in f1bed:
            if b_[0] == c_:
                f1cc[int(b_[1]):int(b_[2])] = 1
                if (int(b_[1]) - last) < filt:
                    f1cc[last:int(b_[1])] = 1
                last = b_[2]

        last = 0
        for b_ in f2bed:
            if b_[0] == c_:
                f2cc[int(b_[1]):int(b_[2])] = 1
                if (int(b_[1]) - last) < filt:
                    f2cc[last:int(b_[1])] = 1
                last = b_[2]

        last = 0
        for b_ in f3bed:
            if b_[0] == c_:
                f3cc[int(b_[1]):int(b_[2])] = 1
                if (int(b_[1]) - last) < filt:
                    f3cc[last:int(b_[1])] = 1
                last = b_[2]

        last = 0
        for b_ in f4bed:
            if b_[0] == c_:
                f4cc[int(b_[1]):int(b_[2])] = 1
                if (int(b_[1]) - last) < filt:
                    f4cc[last:int(b_[1])] = 1
                last = b_[2]

        print(c_)
        cc_correl[c_] = [correl(f1cc, f1cc), correl(f1cc, f2cc), correl(f1cc, f3cc), correl(f1cc, f4cc)]

    return cc_correl


f1name = '/Users/mgabitto/Dropbox (Simons Foundation)/2_PeakCalling/ATACseq_Data/10x_scTest/GM12878_mix_mm10/hg19/atac_v1_hg_10k_peaks.bed'
f2name = '/Users/mgabitto/Dropbox (Simons Foundation)/2_PeakCalling/ATACseq_Data/10x_scTest/GM12878_mix_mm10/hg19-5k/atac_v1_hg_5k_peaks.bed'
f3name = '/Users/mgabitto/Dropbox (Simons Foundation)/2_PeakCalling/ATACseq_Data/10x_scTest/GM12878_mix_mm10/hg19-1k/atac_v1_hg_1k_peaks.bed'
f4name = '/Users/mgabitto/Dropbox (Simons Foundation)/2_PeakCalling/ATACseq_Data/10x_scTest/GM12878_mix_mm10/hg19-500/atac_v1_hg_500_peaks.bed'
bed10x_hgcorrel = bed_file_correlation(f1name, f2name, f3name, f4name, human_chromlens)
pickle.dump(bed10x_hgcorrel, open('human_bed10x_correl.pkl', 'wb'))

f1name = '/Users/mgabitto/Dropbox (Simons Foundation)/2_PeakCalling/ATACseq_Data/10x_scTest/GM12878_mix_mm10/mm10/atac_v1_mm_10k_peaks.bed'
f2name = '/Users/mgabitto/Dropbox (Simons Foundation)/2_PeakCalling/ATACseq_Data/10x_scTest/GM12878_mix_mm10/mm10-5k/atac_v1_mm_5k_peaks.bed'
f3name = '/Users/mgabitto/Dropbox (Simons Foundation)/2_PeakCalling/ATACseq_Data/10x_scTest/GM12878_mix_mm10/mm10-1k/atac_v1_mm_1k_peaks.bed'
f4name = '/Users/mgabitto/Dropbox (Simons Foundation)/2_PeakCalling/ATACseq_Data/10x_scTest/GM12878_mix_mm10/mm10-500/atac_v1_mm_500_peaks.bed'
bed10x_mmcorrel = bed_file_correlation(f1name, f2name, f3name, f4name, mouse_chromlens)
pickle.dump(bed10x_mmcorrel, open('mouse_bed10x_correl.pkl', 'wb'))

hgcor = []
for k_ in bed10x_hgcorrel:
    if not np.isnan(bed10x_hgcorrel[k_][0]):
        hgcor.append(bed10x_hgcorrel[k_])
hgcor.remove(hgcor[14])

mmcor = []
for k_ in bed10x_mmcorrel:
    if not np.isnan(bed10x_mmcorrel[k_][0]):
        mmcor.append(bed10x_mmcorrel[k_])
mmcor.remove(mmcor[9])

fig1 = plt.figure()
plt.plot(np.array(hgcor)[:, :, 0])
pp = PdfPages('10x_hgcor.pdf')
pp.savefig(fig1)
pp.close()

fig1 = plt.figure()
plt.plot(np.array(mmcor)[:, :, 0])
pp = PdfPages('10x_mmcor.pdf')
pp.savefig(fig1)
pp.close()

print(np.mean(mmcor, axis=0))
print(np.mean(hgcor, axis=0))
# ###############################################################
# ###############################################################
# ###############################################################

# ###############################################################
# ###############################################################
# ###############################################################


# COMPUTE 855 BED CORRELATION
f1name = '/Users/mgabitto/Desktop/th17/gb_10k_53b_allpeakssls500.bed'
f2name = '/Users/mgabitto/Desktop/th17/gb_5k_53b_allpeakssls500.bed'
f3name = '/Users/mgabitto/Desktop/th17/gb_1k_53b_allpeakssls500.bed'
f4name = '/Users/mgabitto/Desktop/th17/gb_500_53b_allpeakssls500.bed'
bed53b_mmcorrel = bed_file_correlation(f1name, f2name, f3name, f4name, mouse_chromlens)


f1name = '/Users/mgabitto/Desktop/Projects/PeakCalling/ChromA/mm10-0chr19-10k_allpeaks.bed'
f2name = '/Users/mgabitto/Desktop/Projects/PeakCalling/ChromA/mm10-0chr19-5k_allpeaks.bed'
f3name = '/Users/mgabitto/Desktop/Projects/PeakCalling/ChromA/mm10-0chr19-1k_allpeaks.bed'
f4name = '/Users/mgabitto/Desktop/Projects/PeakCalling/ChromA/mm10-0chr19-500_allpeaks.bed'
bed855_mmcorrel = bed_file_correlation(f1name, f2name, f3name, f4name, mouse_chromlens)

f1name = '/Users/mgabitto/Desktop/th17/mm_10k_53_allpeaks.bed'
f1name = '/Users/mgabitto/Dropbox (Simons Foundation)/2_PeakCalling/ATACseq_Data/10x_scTest/GM12878_mix_mm10/mm10/atac_v1_mm_10k_peaks.bed'
f1name = '/Users/mgabitto/Dropbox (Simons Foundation)/2_PeakCalling/ATACseq_Data/10x_scTest/GM12878_mix_mm10/mm10-500/atac_v1_mm_500_peaks.bed'
f2name = '/Users/mgabitto/Desktop/th17/mm_5k_53_allpeaks.bed'
f3name = '/Users/mgabitto/Desktop/th17/mm_1k_53_allpeaks.bed'
f4name = '/Users/mgabitto/Desktop/th17/mm_500_53_allpeaks.bed'
bed4 = bed_file_correlation(f1name, f2name, f3name, f4name, mouse_chromlens)


f1name = '/Users/mgabitto/Desktop/th17/gb_10k_53_allpeaks.bed'
f2name = '/Users/mgabitto/Desktop/th17/gb_5k_53_allpeaks.bed'
f3name = '/Users/mgabitto/Desktop/th17/gb_1k_53_allpeaks.bed'
f4name = '/Users/mgabitto/Desktop/th17/gb_500_53_allpeaks.bed'
bed53_hgcorrelB = bed_file_correlation(f1name, f2name, f3name, f4name, human_chromlens, filt=500)


f1name = '/Users/mgabitto/Desktop/PD/hg19_10k_rounded.pd'
f2name = '/Users/mgabitto/Desktop/PD/hg19_5k_rounded.pd'
f3name = '/Users/mgabitto/Desktop/PD/hg19_1k_rounded.pd'
f4name = '/Users/mgabitto/Desktop/PD/hg19_500_rounded.pd'
bedPD_hgcorrel = bed_file_correlation(f1name, f2name, f3name, f4name, human_chromlens)


# ###############################################################
# ###############################################################
# ###############################################################

# ###############################################################
# ###############################################################
# ###############################################################
# COMPUTE CLASS TYPE
def get_overlap(a, b):
    return max(0, min(a[1], b[1]) - max(a[0], b[0]))


def peak_class(name_bedh, name_bedl, name_prom):
    bedh = dh.read_bed(name_bedh)
    bedl = dh.read_bed(name_bedl)
    prom = dh.read_bed(name_prom)

    count2 = 0
    count1 = 0
    count0 = 0
    for i_, p_ in enumerate(prom):
        flag = 0
        print(i_ / len(prom))
        for test in bedh:
            if test[0] == p_[0]:
                if getOverlap(test[1:3], p_[1:3]) > 0:
                    count2 += 1
                    flag = 1
                    break

        if flag == 0:
            for test in bedl:
                if test[0] == p_[0]:
                    if getOverlap(test[1:3], p_[1:3]) > 0:
                        count1 += 1
                        flag = 1
                        break
        if flag == 0:
            count0 += 1

    assert(count0 + count1 + count2 == len(prom))

    return [count0, count1, count2]


nhigh = '/Users/mgabitto/Desktop/temp/scmm/gb_10k_1055_2_0.05.bed'
nlow = '/Users/mgabitto/Desktop/temp/scmm/gb_10k_1055_1_0.05.bed'
nprom = '/Users/mgabitto/Desktop/hg19_prom_5k_1k.bed'
count_class_human = peak_class(nhigh, nlow, nprom)

nhigh = '/Users/mgabitto/Desktop/temp/scmm/mm_10k_1055_2_0.05.bed'
nlow = '/Users/mgabitto/Desktop/temp/scmm/mm_10k_1055_1_0.05.bed'
nprom = '/Users/mgabitto/Desktop/Projects/PeakCalling/PeakCalling/Chrom_light2/data/promoters/prom_mm10_genes.bed'
count_class_mouse = peak_class(nhigh, nlow, nprom)
# ###############################################################
# ###############################################################
# ###############################################################


# ###############################################################
# ###############################################################
# ###############################################################
# GENERATE ROBUSTNESS DATASETS
f1name = '/Users/mgabitto/Desktop/MethodsPaper/10x_scTest/GM12878_mix_mm10/hg19/atac_v1_hg_10k_peaks.bed'
hg_chr19 = dh.chr_reads(f1name, 'chr22', 1, human_chromlens['chr22'])
hg_prom_chr19 = dh.read_bed('data/promoters/prom_hg19_genes.bed')
hg_prom_chr19vector = np.zeros(human_chromlens['chr19'])
for coord in hg_prom_chr19:
    if coord[0] == 'chr19':
        hg_prom_chr19vector[b[1]:b[2]] = 1

# Corrupt Read Depth


# Corrupt TSS Enrichment
poss_idx = np.where(hg_prom_chr19vector == 0)[0]







# ###############################################################
# ###############################################################
# ###############################################################
def chip_seq_genome_states(chipseq, peaks):
    from pybedtools import BedTool as bt

    cs = bt(chipseq)

    def len_filter(feature, chr):
        return (feature.chrom == 'chr' + chr)

    frac = []
    length = []
    count = []

    for p_ in peaks:
        bd = bt(p_)
        # bd_and_cs = bd.intersect(cs)
        # count.append(len(bd_and_cs))
        temp = 0
        for i in np.arange(21, 22):
            temp += cs.filter(len_filter, i.__str__()).intersect(bd).count()
        count.append(temp)
        length.append(len(bd))
        print(count)

        cc = {}
        for key in chromlens:
            cc[key] = 0
        for interval in bd:
            cc[interval[0]] += np.abs(np.int(interval[2]) - np.int(interval[1]))

        temp = []
        for key in chromlens:
            temp.append(np.float(cc[key]) / chromlens[key])
        frac.append(np.mean(temp))

    return count, length, frac


def plot_cs(wt_count, wt_length, wt_frac):
    fig401 = plt.figure()
    ax = fig401.subplots(3, 1)
    ax[0].bar(['10x', 'l2', 'l3', 'l4', 'h5', 'h6', 'h7', 'h8', 'h9', 'h10', 'h11', 'h12', 'h13', 'h14', 'h15', 'h16'],
              wt_count, width=0.7)
    plt.xlabel('Algorithm')
    plt.ylabel('Chipseq Recalled')

    ax[1].bar(['10x', 'l2', 'l3', 'l4', 'h5', 'h6', 'h7', 'h8', 'h9', 'h10', 'h11', 'h12', 'h13', 'h14', 'h15', 'h16'],
              wt_length, width=0.7)
    plt.xlabel('Algorithm')
    plt.ylabel('Number of peaks')

    ax[2].bar(['10x', 'l2', 'l3', 'l4', 'h5', 'h6', 'h7', 'h8', 'h9', 'h10', 'h11', 'h12', 'h13', 'h14', 'h15', 'h16'],
              wt_frac, width=0.7)
    plt.xlabel('Algorithm')
    plt.ylabel('Accessible Genome Fraction')


# NUMBER OF CHIP-SEQ PEAKS CAPTURED
chromlens = {'chr1': 249250621, 'chr2': 243199373, 'chr3': 198022430, 'chr4': 191154276, 'chr5': 180915260,
              'chr6': 171115067, 'chr7': 159138663, 'chr8': 146364022, 'chr9': 141213431, 'chr10': 135534747,
              'chr11': 135006516, 'chr12': 133851895, 'chr13': 115169878, 'chr14': 107349540,
              'chr15': 102531392, 'chr16': 90354753, 'chr17': 81195210, 'chr18': 78077248, 'chr19': 59128983,
              'chr20': 63025520, 'chr22': 51304566, 'chr21': 48129895, 'chrX': 155270560, 'chrY': 59373566,
              'chrM': 16571}

from pybedtools import BedTool as bt
# chipseq = '/Users/mgabitto/Dropbox (Simons Foundation)/2_PeakCalling/Chipseq/ReMapTfsGM12878/remap2018_GM12878_summits.bed'
chipseq = '/Users/mgabitto/Dropbox (Simons Foundation)/2_PeakCalling/Chipseq/ReMapTfsGM12878/b.bed'

ll = []
ll.append('/Users/mgabitto/Desktop/th17/gb_10k_53_allpeaks.bed')
ll.append('/Users/mgabitto/Desktop/th17/gb_5k_53_allpeaks.bed')
ll.append('/Users/mgabitto/Desktop/th17/gb_1k_53_allpeaks.bed')
ll.append('/Users/mgabitto/Desktop/th17/gb_500_53_allpeaks.bed')
ll.append('/Users/mgabitto/Desktop/PD/hg19_10k_roundedtab.pd')
ll.append('/Users/mgabitto/Desktop/PD/hg19_5k_roundedtab.pd')
ll.append('/Users/mgabitto/Desktop/PD/hg19_1k_roundedtab.pd')
ll.append('/Users/mgabitto/Desktop/PD/hg19_500_roundedtab.pd')
ll.append('/Users/mgabitto/Dropbox (Simons Foundation)/2_PeakCalling/ATACseq_Data/10x_scTest/GM12878_mix_mm10/hg19/atac_v1_hg_10k_peaks.bed')
ll.append('/Users/mgabitto/Dropbox (Simons Foundation)/2_PeakCalling/ATACseq_Data/10x_scTest/GM12878_mix_mm10/hg19-5k/atac_v1_hg_5k_peaks.bed')
ll.append('/Users/mgabitto/Dropbox (Simons Foundation)/2_PeakCalling/ATACseq_Data/10x_scTest/GM12878_mix_mm10/hg19-1k/atac_v1_hg_1k_peaks.bed')
ll.append('/Users/mgabitto/Dropbox (Simons Foundation)/2_PeakCalling/ATACseq_Data/10x_scTest/GM12878_mix_mm10/hg19-500/atac_v1_hg_500_peaks.bed')
ll.append('/Users/mgabitto/Desktop/th17/gb_10k_75_allpeaks.bed')
ll.append('/Users/mgabitto/Desktop/th17/gb_5k_75_allpeaks.bed')
ll.append('/Users/mgabitto/Desktop/th17/gb_1k_75_allpeaks.bed')
ll.append('/Users/mgabitto/Desktop/th17/gb_500_75_allpeaks.bed')
ll.append('/Users/mgabitto/Desktop/temp/MACS/peaks_hg19f.bed')


ll.append('/Users/mgabitto/Desktop/th17/gb_10k_75_allpeakss500.bed')
ll.append('/Users/mgabitto/Desktop/th17/gb_5k_75_allpeakss500.bed')
ll.append('/Users/mgabitto/Desktop/th17/gb_1k_75_allpeakss500.bed')
ll.append('/Users/mgabitto/Desktop/th17/gb_500_75_allpeakss500.bed')
ll.append('/Users/mgabitto/Desktop/th17/gb_10k_53_allpeakss500.bed')
ll.append('/Users/mgabitto/Desktop/th17/gb_5k_53_allpeakss500.bed')
ll.append('/Users/mgabitto/Desktop/th17/gb_1k_53_allpeakss500.bed')
ll.append('/Users/mgabitto/Desktop/th17/gb_500_53_allpeakss500.bed')

ll = []
ll.append('/Users/mgabitto/Desktop/th17/gb_10k_75b_allpeakss500.bed')
ll.append('/Users/mgabitto/Desktop/th17/gb_5k_75b_allpeakss500.bed')
ll.append('/Users/mgabitto/Desktop/th17/gb_1k_75b_allpeakss500.bed')
ll.append('/Users/mgabitto/Desktop/th17/gb_500_75b_allpeakss500.bed')
ll.append('/Users/mgabitto/Desktop/th17/gb_10k_53b_allpeakss500.bed')
ll.append('/Users/mgabitto/Desktop/th17/gb_5k_53b_allpeakss500.bed')
ll.append('/Users/mgabitto/Desktop/th17/gb_1k_53b_allpeakss500.bed')
ll.append('/Users/mgabitto/Desktop/th17/gb_500_53b_allpeakss500.bed')

w_count, w_length, w_frac = chip_seq_genome_states(chipseq, ll)
w_count2, w_length2, w_frac2 = chip_seq_genome_states(chipseq, ll2)
plot_cs(w_count, w_length, w_frac)

ll = []
ll.append('/Users/mgabitto/Desktop/th17/gb_10k_53b_allpeakssls500.bed')
ll.append('/Users/mgabitto/Desktop/th17/gb_5k_53b_allpeakssls500.bed')
ll.append('/Users/mgabitto/Desktop/th17/gb_1k_53b_allpeakssls500.bed')
ll.append('/Users/mgabitto/Desktop/th17/gb_500_53b_allpeakssls500.bed')
w_count21, w_length21, w_frac21 = chip_seq_genome_states(chipseq, ll)

chipseq = '/Users/mgabitto/Dropbox (Simons Foundation)/2_PeakCalling/Chipseq/ReMapTfsGM12878/b.bed'
bedcs = dh.read_bed(chipseq)

test ='/Users/mgabitto/Desktop/th17/gb_500_53b_allpeakssls500.bed'
test = '/Users/mgabitto/Dropbox (Simons Foundation)/2_PeakCalling/ATACseq_Data/10x_scTest/GM12878_mix_mm10/hg19-500/atac_v1_hg_500_peaks.bed'

bedtest = dh.read_bed(test)
bedtesta = np.array(bedtest)
hist = np.int32(bedtesta[:, 2]) - np.int32(bedtesta[:, 1])

fig1 = plt.figure()
plt.hist(hist, 100)
pp = PdfPages('hist50010x.pdf')
pp.savefig(fig1)
pp.close()

