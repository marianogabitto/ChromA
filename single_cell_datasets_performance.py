import classes.data_handle as dh
import numpy as np
import pickle
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages


# #####################################################################
# #####################################################################
# #####################################################################
# COMPUTE CORRELATION
def correl(a, b):
    a1 = (a - np.mean(a)) / (np.std(a) * len(a))
    b1 = (b - np.mean(b)) / (np.std(b))

    return np.correlate(a1[:, 0], b1[:, 0])


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


human_chromlens = dh.human_lens()
human_chromlens.__delitem__('chrM')
human_chromlens.__delitem__('chrX')
human_chromlens.__delitem__('chrY')
mouse_chromlens = dh.mouse_lens()
mouse_chromlens.__delitem__('chrM')
mouse_chromlens.__delitem__('chrX')
mouse_chromlens.__delitem__('chrY')


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

f1name = '/Users/mgabitto/Desktop/th17/gb_10k_53b_allpeakssls500.bed'
f2name = '/Users/mgabitto/Desktop/th17/gb_5k_53b_allpeakssls500.bed'
f3name = '/Users/mgabitto/Desktop/th17/gb_1k_53b_allpeakssls500.bed'
f4name = '/Users/mgabitto/Desktop/th17/gb_500_53b_allpeakssls500.bed'
bed53_hgcorrel = bed_file_correlation(f1name, f2name, f3name, f4name, human_chromlens)
pickle.dump(bed53_hgcorrel, open('human_bed_correl.pkl', 'wb'))

f1name = '/Users/mgabitto/Desktop/th17/gb_10k_75b_allpeakssls500.bed'
f2name = '/Users/mgabitto/Desktop/th17/gb_5k_75b_allpeakssls500.bed'
f3name = '/Users/mgabitto/Desktop/th17/gb_1k_75b_allpeakssls500.bed'
f4name = '/Users/mgabitto/Desktop/th17/gb_500_75b_allpeakssls500.bed'
bed75_hgcorrel = bed_file_correlation(f1name, f2name, f3name, f4name, human_chromlens)
pickle.dump(bed75_hgcorrel, open('human_bed_correl.pkl', 'wb'))

f1name = '/Users/mgabitto/Desktop/th17/mm_10k_53b_allpeakssls500.bed'
f2name = '/Users/mgabitto/Desktop/th17/mm_5k_53b_allpeakssls500.bed'
f3name = '/Users/mgabitto/Desktop/th17/mm_1k_53b_allpeakssls500.bed'
f4name = '/Users/mgabitto/Desktop/th17/mm_500_53b_allpeakssls500.bed'
bed53_mmcorrel = bed_file_correlation(f1name, f2name, f3name, f4name, mouse_chromlens)
pickle.dump(bed53_mmcorrel, open('mouse_bed_correl.pkl', 'wb'))



f1name = '/Users/mgabitto/Desktop/PD/hg19_10k_roundedtab.pd'
f2name = '/Users/mgabitto/Desktop/PD/hg19_5k_roundedtab.pd'
f3name = '/Users/mgabitto/Desktop/PD/hg19_1k_roundedtab.pd'
f4name = '/Users/mgabitto/Desktop/PD/hg19_500_roundedtab.pd'
bedPD_hgcorrel = bed_file_correlation(f1name, f2name, f3name, f4name, human_chromlens)
pickle.dump(bedPD_hgcorrel, open('human_PDbed_correl.pkl', 'wb'))

f1name = '/Users/mgabitto/Desktop/PD/mm10_10k_roundedtab.pd'
f2name = '/Users/mgabitto/Desktop/PD/mm10_5k_roundedtab.pd'
f3name = '/Users/mgabitto/Desktop/PD/mm10_1k_roundedtab.pd'
f4name = '/Users/mgabitto/Desktop/PD/mm10_500_roundedtab.pd'
bedPD_mmcorrel = bed_file_correlation(f1name, f2name, f3name, f4name, mouse_chromlens)
pickle.dump(bedPD_mmcorrel, open('mouse_PDbed_correl.pkl', 'wb'))

f1name = '/Users/mgabitto/Desktop/MACS/10k/atac_v1_hg_10k_possorted_bam_peaks.bed'
f2name = '/Users/mgabitto/Desktop/MACS/5k/atac_v1_hg_5k_possorted_bam_peaks.bed'
f3name = '/Users/mgabitto/Desktop/MACS/1k/atac_v1_hg_1k_possorted_bam_peaks.bed'
f4name = '/Users/mgabitto/Desktop/MACS/500/atac_v1_hg_500_possorted_bam_peaks.bed'
bedM_hgcorrel = bed_file_correlation(f1name, f2name, f3name, f4name, human_chromlens)
pickle.dump(bedM_hgcorrel, open('results/human_Mbed_correl.pkl', 'wb'))

f1name = '/Users/mgabitto/Desktop/MACS/10k/atac_v1_mm_10k_possorted_bam_peaks.bed'
f2name = '/Users/mgabitto/Desktop/MACS/5k/atac_v1_mm_5k_possorted_bam_peaks.bed'
f3name = '/Users/mgabitto/Desktop/MACS/1k/atac_v1_mm_1k_possorted_bam_peaks.bed'
f4name = '/Users/mgabitto/Desktop/MACS/500/atac_v1_mm_500_possorted_bam_peaks.bed'
bedM_mmcorrel = bed_file_correlation(f1name, f2name, f3name, f4name, mouse_chromlens)
pickle.dump(bedM_mmcorrel, open('results/mouse_Mbed_correl.pkl', 'wb'))


fig1 = plt.figure()
plt.plot(np.array(hgcor)[:, :, 0])
pp = PdfPages('10x_hgcor.pdf')
pp.savefig(fig1)
pp.close()

bed10x_hgcorrel = pickle.load(open('/Users/mgabitto/Desktop/Projects/PeakCalling/ChromA/results/human_bed10x_correl.pkl', 'rb'))
res10xhg = np.zeros((1, 4))
hgcor = []
for k_ in bed10x_hgcorrel.keys():
    if not np.any(np.isnan(bed10x_hgcorrel[k_])):
        hgcor.append(bed10x_hgcorrel[k_])
res10xhg[0, :] = np.mean(np.array(hgcor), axis=0)[:, 0]


# ###############################################################
# ###############################################################
# ###############################################################


# ###############################################################
# ###############################################################
# ###############################################################
from pybedtools import BedTool as bt


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
        for i in np.arange(1, 22):
            temp += cs.filter(len_filter, i.__str__()).intersect(bd).count()
        count.append(temp)
        length.append(len(bd))
        print(count)

        """
        cc = {}
        for key in chromlens:
            cc[key] = 0
        for interval in bd:
            cc[interval[0]] += np.abs(np.int(interval[2]) - np.int(interval[1]))

        temp = []
        for key in chromlens:
            temp.append(np.float(cc[key]) / chromlens[key])
        frac.append(np.mean(temp))
        """
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
# chipseq = '/Users/mgabitto/Dropbox (Simons Foundation)/2_PeakCalling/Chipseq/ReMapTfsGM12878/remap2018_GM12878_summits.bed'
chipseq = '/Users/mgabitto/Dropbox (Simons Foundation)/2_PeakCalling/Chipseq/ReMapTfsGM12878/b.bed'

ll = []
ll.append('/Users/mgabitto/Dropbox (Simons Foundation)/2_PeakCalling/ATACseq_Data/10x_scTest/GM12878_mix_mm10/hg19/atac_v1_hg_10k_peaks.bed')
ll.append('/Users/mgabitto/Dropbox (Simons Foundation)/2_PeakCalling/ATACseq_Data/10x_scTest/GM12878_mix_mm10/hg19-5k/atac_v1_hg_5k_peaks.bed')
ll.append('/Users/mgabitto/Dropbox (Simons Foundation)/2_PeakCalling/ATACseq_Data/10x_scTest/GM12878_mix_mm10/hg19-1k/atac_v1_hg_1k_peaks.bed')
ll.append('/Users/mgabitto/Dropbox (Simons Foundation)/2_PeakCalling/ATACseq_Data/10x_scTest/GM12878_mix_mm10/hg19-500/atac_v1_hg_500_peaks.bed')
ll.append('/Users/mgabitto/Desktop/PD/hg19_10k_roundedtab.pd')
ll.append('/Users/mgabitto/Desktop/PD/hg19_5k_roundedtab.pd')
ll.append('/Users/mgabitto/Desktop/PD/hg19_1k_roundedtab.pd')
ll.append('/Users/mgabitto/Desktop/PD/hg19_500_roundedtab.pd')
ll.append('/Users/mgabitto/Desktop/MACS/10k/atac_v1_hg_10k_possorted_bam_peaks.bed')
ll.append('/Users/mgabitto/Desktop/MACS/5k/atac_v1_hg_5k_possorted_bam_peaks.bed')
ll.append('/Users/mgabitto/Desktop/MACS/1k/atac_v1_hg_1k_possorted_bam_peaks.bed')
ll.append('/Users/mgabitto/Desktop/MACS/500/atac_v1_hg_500_possorted_bam_peaks.bed')

chromlens = dh.human_lens()
w_count, w_length, w_frac = chip_seq_genome_states(chipseq, ll)

count_chr1_21_10x = w_count[0:4]
count_chr1_21_PD = w_count[4:8]
count_chr1_21_M = w_count[8:12]
np.save('results/count_chr1_21_10x.npy', count_chr1_21_10x)
np.save('results/count_chr1_21_PD.npy', count_chr1_21_PD)
np.save('results/count_chr1_21_M.npy', count_chr1_21_M)
# ###############################################################
# ###############################################################
# ###############################################################

# ###############################################################
# ###############################################################
# ###############################################################
chipseq = '/Users/mgabitto/Dropbox (Simons Foundation)/2_PeakCalling/Chipseq/ReMapTfsGM12878/b.bed'

ll = []
ll.append('/Users/mgabitto/Desktop/th17/gb_sub/gb_10k_53b_allpeaks_100sl_250s.bed')
ll.append('/Users/mgabitto/Desktop/th17/gb_sub/gb_5k_53b_allpeaks_100sl_250s.bed')
ll.append('/Users/mgabitto/Desktop/th17/gb_sub/gb_1k_53b_allpeaks_100sl_250s.bed')
ll.append('/Users/mgabitto/Desktop/th17/gb_sub/gb_500_53b_allpeaks_100sl_250s.bed')

chromlens = dh.human_lens()
w_count_53_100_250, w_length, w_frac = chip_seq_genome_states(chipseq, ll)
np.save('results/count_chr1_21_53_100_250.npy', w_count_53_100_250)


ll = []
ll.append('/Users/mgabitto/Desktop/th17/gb_sub/gb_10k_53b_allpeaks_100sl_500s.bed')
ll.append('/Users/mgabitto/Desktop/th17/gb_sub/gb_5k_53b_allpeaks_100sl_500s.bed')
ll.append('/Users/mgabitto/Desktop/th17/gb_sub/gb_1k_53b_allpeaks_100sl_500s.bed')
ll.append('/Users/mgabitto/Desktop/th17/gb_sub/gb_500_53b_allpeaks_100sl_500s.bed')

chromlens = dh.human_lens()
w_count_53_100_500, w_length, w_frac = chip_seq_genome_states(chipseq, ll)
np.save('results/count_chr1_21_53_100_500.npy', w_count_53_100_500)








ll = []
ll.append('/Users/mgabitto/Desktop/th17/gb_sub/gb_10k_53b_allpeaks_50sl_250s.bed')
ll.append('/Users/mgabitto/Desktop/th17/gb_sub/gb_5k_53b_allpeaks_50sl_250s.bed')
ll.append('/Users/mgabitto/Desktop/th17/gb_sub/gb_1k_53b_allpeaks_50sl_250s.bed')
ll.append('/Users/mgabitto/Desktop/th17/gb_sub/gb_500_53b_allpeaks_50sl_250s.bed')

chromlens = dh.human_lens()
w_count_53_50_250, w_length, w_frac = chip_seq_genome_states(chipseq, ll)
np.save('results/count_chr1_21_53_50_250.npy', w_count_53_50_250)


ll = []
ll.append('/Users/mgabitto/Desktop/th17/gb_sub/gb_10k_53b_allpeaks_50sl_500s.bed')
ll.append('/Users/mgabitto/Desktop/th17/gb_sub/gb_5k_53b_allpeaks_50sl_500s.bed')
ll.append('/Users/mgabitto/Desktop/th17/gb_sub/gb_1k_53b_allpeaks_50sl_500s.bed')
ll.append('/Users/mgabitto/Desktop/th17/gb_sub/gb_500_53b_allpeaks_50sl_500s.bed')

chromlens = dh.human_lens()
w_count_53_50_500, w_length, w_frac = chip_seq_genome_states(chipseq, ll)
np.save('results/count_chr1_21_53_50_500.npy', w_count_53_50_500)







ll = []
ll.append('/Users/mgabitto/Desktop/th17/gb_sub/gb_10k_53b_allpeaks_25sl_250s.bed')
ll.append('/Users/mgabitto/Desktop/th17/gb_sub/gb_5k_53b_allpeaks_25sl_250s.bed')
ll.append('/Users/mgabitto/Desktop/th17/gb_sub/gb_1k_53b_allpeaks_25sl_250s.bed')
ll.append('/Users/mgabitto/Desktop/th17/gb_sub/gb_500_53b_allpeaks_25sl_250s.bed')

chromlens = dh.human_lens()
w_count_53_25_250, w_length, w_frac = chip_seq_genome_states(chipseq, ll)
np.save('results/count_chr1_21_53_50_250.npy', w_count_53_25_250)


ll = []
ll.append('/Users/mgabitto/Desktop/th17/gb_sub/gb_10k_53b_allpeaks_25sl_500s.bed')
ll.append('/Users/mgabitto/Desktop/th17/gb_sub/gb_5k_53b_allpeaks_25sl_500s.bed')
ll.append('/Users/mgabitto/Desktop/th17/gb_sub/gb_1k_53b_allpeaks_25sl_500s.bed')
ll.append('/Users/mgabitto/Desktop/th17/gb_sub/gb_500_53b_allpeaks_25sl_500s.bed')

chromlens = dh.human_lens()
w_count_53_25_500, w_length, w_frac = chip_seq_genome_states(chipseq, ll)
np.save('results/count_chr1_21_53_25_500.npy', w_count_53_25_500)




# ###############################################################
# ###############################################################
# ###############################################################
# ###############################################################
# ###############################################################
# ###############################################################
chipseq = '/Users/mgabitto/Dropbox (Simons Foundation)/2_PeakCalling/Chipseq/ReMapTfsGM12878/b.bed'

ll = []
ll.append('/Users/mgabitto/Desktop/th17/gb_sub/gb_10k_75b_allpeaks_100sl_250s.bed')
ll.append('/Users/mgabitto/Desktop/th17/gb_sub/gb_5k_75b_allpeaks_100sl_250s.bed')
ll.append('/Users/mgabitto/Desktop/th17/gb_sub/gb_1k_75b_allpeaks_100sl_250s.bed')
ll.append('/Users/mgabitto/Desktop/th17/gb_sub/gb_500_75b_allpeaks_100sl_250s.bed')

chromlens = dh.human_lens()
w_count_75_100_250, w_length, w_frac = chip_seq_genome_states(chipseq, ll)
np.save('results/count_chr1_21_75_100_250.npy', w_count_75_100_250)


ll = []
ll.append('/Users/mgabitto/Desktop/th17/gb_sub/gb_10k_75b_allpeaks_100sl_500s.bed')
ll.append('/Users/mgabitto/Desktop/th17/gb_sub/gb_5k_75b_allpeaks_100sl_500s.bed')
ll.append('/Users/mgabitto/Desktop/th17/gb_sub/gb_1k_75b_allpeaks_100sl_500s.bed')
ll.append('/Users/mgabitto/Desktop/th17/gb_sub/gb_500_75b_allpeaks_100sl_500s.bed')

chromlens = dh.human_lens()
w_count_75_100_500, w_length, w_frac = chip_seq_genome_states(chipseq, ll)
np.save('results/count_chr1_21_75_100_500.npy', w_count_75_100_500)








ll = []
ll.append('/Users/mgabitto/Desktop/th17/gb_sub/gb_10k_75b_allpeaks_50sl_250s.bed')
ll.append('/Users/mgabitto/Desktop/th17/gb_sub/gb_5k_75b_allpeaks_50sl_250s.bed')
ll.append('/Users/mgabitto/Desktop/th17/gb_sub/gb_1k_75b_allpeaks_50sl_250s.bed')
ll.append('/Users/mgabitto/Desktop/th17/gb_sub/gb_500_75b_allpeaks_50sl_250s.bed')

chromlens = dh.human_lens()
w_count_75_50_250, w_length, w_frac = chip_seq_genome_states(chipseq, ll)
np.save('results/count_chr1_21_75_50_250.npy', w_count_75_50_250)


ll = []
ll.append('/Users/mgabitto/Desktop/th17/gb_sub/gb_10k_75b_allpeaks_50sl_500s.bed')
ll.append('/Users/mgabitto/Desktop/th17/gb_sub/gb_5k_75b_allpeaks_50sl_500s.bed')
ll.append('/Users/mgabitto/Desktop/th17/gb_sub/gb_1k_75b_allpeaks_50sl_500s.bed')
ll.append('/Users/mgabitto/Desktop/th17/gb_sub/gb_500_75b_allpeaks_50sl_500s.bed')

chromlens = dh.human_lens()
w_count_75_50_500, w_length, w_frac = chip_seq_genome_states(chipseq, ll)
np.save('results/count_chr1_21_75_50_500.npy', w_count_75_50_500)







ll = []
ll.append('/Users/mgabitto/Desktop/th17/gb_sub/gb_10k_75b_allpeaks_25sl_250s.bed')
ll.append('/Users/mgabitto/Desktop/th17/gb_sub/gb_5k_75b_allpeaks_25sl_250s.bed')
ll.append('/Users/mgabitto/Desktop/th17/gb_sub/gb_1k_75b_allpeaks_25sl_250s.bed')
ll.append('/Users/mgabitto/Desktop/th17/gb_sub/gb_500_75b_allpeaks_25sl_250s.bed')

chromlens = dh.human_lens()
w_count_75_25_250, w_length, w_frac = chip_seq_genome_states(chipseq, ll)
np.save('results/count_chr1_21_75_50_250.npy', w_count_75_25_250)


ll = []
ll.append('/Users/mgabitto/Desktop/th17/gb_sub/gb_10k_75b_allpeaks_25sl_500s.bed')
ll.append('/Users/mgabitto/Desktop/th17/gb_sub/gb_5k_75b_allpeaks_25sl_500s.bed')
ll.append('/Users/mgabitto/Desktop/th17/gb_sub/gb_1k_75b_allpeaks_25sl_500s.bed')
ll.append('/Users/mgabitto/Desktop/th17/gb_sub/gb_500_75b_allpeaks_25sl_500s.bed')

chromlens = dh.human_lens()
w_count_75_25_500, w_length, w_frac = chip_seq_genome_states(chipseq, ll)
np.save('results/count_chr1_21_75_25_500.npy', w_count_75_25_500)

















# #####################################################
bed53 = []
f1name = '/Users/mgabitto/Desktop/th17/gb_sub/gb_10k_53b_allpeaks_100sl_250s.bed'
f2name = '/Users/mgabitto/Desktop/th17/gb_sub/gb_5k_53b_allpeaks_100sl_250s.bed'
f3name = '/Users/mgabitto/Desktop/th17/gb_sub/gb_1k_53b_allpeaks_100sl_250s.bed'
f4name = '/Users/mgabitto/Desktop/th17/gb_sub/gb_500_53b_allpeaks_100sl_250s.bed'
bed53.append(bed_file_correlation(f1name, f2name, f3name, f4name, human_chromlens))

f1name = '/Users/mgabitto/Desktop/th17/gb_sub/gb_10k_53b_allpeaks_100sl_500s.bed'
f2name = '/Users/mgabitto/Desktop/th17/gb_sub/gb_5k_53b_allpeaks_100sl_500s.bed'
f3name = '/Users/mgabitto/Desktop/th17/gb_sub/gb_1k_53b_allpeaks_100sl_500s.bed'
f4name = '/Users/mgabitto/Desktop/th17/gb_sub/gb_500_53b_allpeaks_100sl_500s.bed'
bed53.append(bed_file_correlation(f1name, f2name, f3name, f4name, human_chromlens))

f1name = '/Users/mgabitto/Desktop/th17/gb_sub/gb_10k_53b_allpeaks_50sl_250s.bed'
f2name = '/Users/mgabitto/Desktop/th17/gb_sub/gb_5k_53b_allpeaks_50sl_250s.bed'
f3name = '/Users/mgabitto/Desktop/th17/gb_sub/gb_1k_53b_allpeaks_50sl_250s.bed'
f4name = '/Users/mgabitto/Desktop/th17/gb_sub/gb_500_53b_allpeaks_50sl_250s.bed'
bed53.append(bed_file_correlation(f1name, f2name, f3name, f4name, human_chromlens))

f1name = '/Users/mgabitto/Desktop/th17/gb_sub/gb_10k_53b_allpeaks_50sl_500s.bed'
f2name = '/Users/mgabitto/Desktop/th17/gb_sub/gb_5k_53b_allpeaks_50sl_500s.bed'
f3name = '/Users/mgabitto/Desktop/th17/gb_sub/gb_1k_53b_allpeaks_50sl_500s.bed'
f4name = '/Users/mgabitto/Desktop/th17/gb_sub/gb_500_53b_allpeaks_50sl_500s.bed'
bed53.append(bed_file_correlation(f1name, f2name, f3name, f4name, human_chromlens))

f1name = '/Users/mgabitto/Desktop/th17/gb_sub/gb_10k_53b_allpeaks_25sl_250s.bed'
f2name = '/Users/mgabitto/Desktop/th17/gb_sub/gb_5k_53b_allpeaks_25sl_250s.bed'
f3name = '/Users/mgabitto/Desktop/th17/gb_sub/gb_1k_53b_allpeaks_25sl_250s.bed'
f4name = '/Users/mgabitto/Desktop/th17/gb_sub/gb_500_53b_allpeaks_25sl_250s.bed'
bed53.append(bed_file_correlation(f1name, f2name, f3name, f4name, human_chromlens))

f1name = '/Users/mgabitto/Desktop/th17/gb_sub/gb_10k_53b_allpeaks_25sl_500s.bed'
f2name = '/Users/mgabitto/Desktop/th17/gb_sub/gb_5k_53b_allpeaks_25sl_500s.bed'
f3name = '/Users/mgabitto/Desktop/th17/gb_sub/gb_1k_53b_allpeaks_25sl_500s.bed'
f4name = '/Users/mgabitto/Desktop/th17/gb_sub/gb_500_53b_allpeaks_25sl_500s.bed'
bed53.append(bed_file_correlation(f1name, f2name, f3name, f4name, human_chromlens))

res53hg = np.zeros((len(bed53), 4))
for i in np.arange(len(bed53)):
    hgcor = []
    for k_ in bed53[i].keys():
        if not np.any(np.isnan(bed53[i][k_])):
            hgcor.append(bed53[i][k_])
    res53hg[i, :] = np.mean(np.array(hgcor), axis=0)[:, 0]
np.save('results/bed53_correl_explorehg.npy', res53hg)



# #####################################################
bed75 = []
f1name = '/Users/mgabitto/Desktop/th17/gb_sub/gb_10k_75b_allpeaks_100sl_250s.bed'
f2name = '/Users/mgabitto/Desktop/th17/gb_sub/gb_5k_75b_allpeaks_100sl_250s.bed'
f3name = '/Users/mgabitto/Desktop/th17/gb_sub/gb_1k_75b_allpeaks_100sl_250s.bed'
f4name = '/Users/mgabitto/Desktop/th17/gb_sub/gb_500_75b_allpeaks_100sl_250s.bed'
bed75.append(bed_file_correlation(f1name, f2name, f3name, f4name, human_chromlens))

f1name = '/Users/mgabitto/Desktop/th17/gb_sub/gb_10k_75b_allpeaks_100sl_500s.bed'
f2name = '/Users/mgabitto/Desktop/th17/gb_sub/gb_5k_75b_allpeaks_100sl_500s.bed'
f3name = '/Users/mgabitto/Desktop/th17/gb_sub/gb_1k_75b_allpeaks_100sl_500s.bed'
f4name = '/Users/mgabitto/Desktop/th17/gb_sub/gb_500_75b_allpeaks_100sl_500s.bed'
bed75.append(bed_file_correlation(f1name, f2name, f3name, f4name, human_chromlens))

f1name = '/Users/mgabitto/Desktop/th17/gb_sub/gb_10k_75b_allpeaks_50sl_250s.bed'
f2name = '/Users/mgabitto/Desktop/th17/gb_sub/gb_5k_75b_allpeaks_50sl_250s.bed'
f3name = '/Users/mgabitto/Desktop/th17/gb_sub/gb_1k_75b_allpeaks_50sl_250s.bed'
f4name = '/Users/mgabitto/Desktop/th17/gb_sub/gb_500_75b_allpeaks_50sl_250s.bed'
bed75.append(bed_file_correlation(f1name, f2name, f3name, f4name, human_chromlens))

f1name = '/Users/mgabitto/Desktop/th17/gb_sub/gb_10k_75b_allpeaks_50sl_500s.bed'
f2name = '/Users/mgabitto/Desktop/th17/gb_sub/gb_5k_75b_allpeaks_50sl_500s.bed'
f3name = '/Users/mgabitto/Desktop/th17/gb_sub/gb_1k_75b_allpeaks_50sl_500s.bed'
f4name = '/Users/mgabitto/Desktop/th17/gb_sub/gb_500_75b_allpeaks_50sl_500s.bed'
bed75.append(bed_file_correlation(f1name, f2name, f3name, f4name, human_chromlens))

f1name = '/Users/mgabitto/Desktop/th17/gb_sub/gb_10k_75b_allpeaks_25sl_250s.bed'
f2name = '/Users/mgabitto/Desktop/th17/gb_sub/gb_5k_75b_allpeaks_25sl_250s.bed'
f3name = '/Users/mgabitto/Desktop/th17/gb_sub/gb_1k_75b_allpeaks_25sl_250s.bed'
f4name = '/Users/mgabitto/Desktop/th17/gb_sub/gb_500_75b_allpeaks_25sl_250s.bed'
bed75.append(bed_file_correlation(f1name, f2name, f3name, f4name, human_chromlens))

f1name = '/Users/mgabitto/Desktop/th17/gb_sub/gb_10k_75b_allpeaks_25sl_500s.bed'
f2name = '/Users/mgabitto/Desktop/th17/gb_sub/gb_5k_75b_allpeaks_25sl_500s.bed'
f3name = '/Users/mgabitto/Desktop/th17/gb_sub/gb_1k_75b_allpeaks_25sl_500s.bed'
f4name = '/Users/mgabitto/Desktop/th17/gb_sub/gb_500_75b_allpeaks_25sl_500s.bed'
bed75.append(bed_file_correlation(f1name, f2name, f3name, f4name, human_chromlens))

res75hg = np.zeros((len(bed75), 4))
for i in np.arange(len(bed75)):
    hgcor = []
    for k_ in bed75[i].keys():
        if not np.any(np.isnan(bed75[i][k_])):
            hgcor.append(bed75[i][k_])
    res75hg[i, :] = np.mean(np.array(hgcor), axis=0)[:, 0]
np.save('results/bed75_correl_explorehg.npy', res75hg)





# #########################################################
bed53mm = []
f1name = '/Users/mgabitto/Desktop/th17/mm_sub/mm_10k_53b_allpeaks_100sl_250s.bed'
f2name = '/Users/mgabitto/Desktop/th17/mm_sub/mm_5k_53b_allpeaks_100sl_250s.bed'
f3name = '/Users/mgabitto/Desktop/th17/mm_sub/mm_1k_53b_allpeaks_100sl_250s.bed'
f4name = '/Users/mgabitto/Desktop/th17/mm_sub/mm_500_53b_allpeaks_100sl_250s.bed'
bed53mm.append(bed_file_correlation(f1name, f2name, f3name, f4name, mouse_chromlens))

f1name = '/Users/mgabitto/Desktop/th17/mm_sub/mm_10k_53b_allpeaks_100sl_500s.bed'
f2name = '/Users/mgabitto/Desktop/th17/mm_sub/mm_5k_53b_allpeaks_100sl_500s.bed'
f3name = '/Users/mgabitto/Desktop/th17/mm_sub/mm_1k_53b_allpeaks_100sl_500s.bed'
f4name = '/Users/mgabitto/Desktop/th17/mm_sub/mm_500_53b_allpeaks_100sl_500s.bed'
bed53mm.append(bed_file_correlation(f1name, f2name, f3name, f4name, mouse_chromlens))

f1name = '/Users/mgabitto/Desktop/th17/mm_sub/mm_10k_53b_allpeaks_50sl_250s.bed'
f2name = '/Users/mgabitto/Desktop/th17/mm_sub/mm_5k_53b_allpeaks_50sl_250s.bed'
f3name = '/Users/mgabitto/Desktop/th17/mm_sub/mm_1k_53b_allpeaks_50sl_250s.bed'
f4name = '/Users/mgabitto/Desktop/th17/mm_sub/mm_500_53b_allpeaks_50sl_250s.bed'
bed53mm.append(bed_file_correlation(f1name, f2name, f3name, f4name, mouse_chromlens))

f1name = '/Users/mgabitto/Desktop/th17/mm_sub/mm_10k_53b_allpeaks_50sl_500s.bed'
f2name = '/Users/mgabitto/Desktop/th17/mm_sub/mm_5k_53b_allpeaks_50sl_500s.bed'
f3name = '/Users/mgabitto/Desktop/th17/mm_sub/mm_1k_53b_allpeaks_50sl_500s.bed'
f4name = '/Users/mgabitto/Desktop/th17/mm_sub/mm_500_53b_allpeaks_50sl_500s.bed'
bed53mm.append(bed_file_correlation(f1name, f2name, f3name, f4name, mouse_chromlens))

f1name = '/Users/mgabitto/Desktop/th17/mm_sub/mm_10k_53b_allpeaks_25sl_250s.bed'
f2name = '/Users/mgabitto/Desktop/th17/mm_sub/mm_5k_53b_allpeaks_25sl_250s.bed'
f3name = '/Users/mgabitto/Desktop/th17/mm_sub/mm_1k_53b_allpeaks_25sl_250s.bed'
f4name = '/Users/mgabitto/Desktop/th17/mm_sub/mm_500_53b_allpeaks_25sl_250s.bed'
bed53mm.append(bed_file_correlation(f1name, f2name, f3name, f4name, mouse_chromlens))

f1name = '/Users/mgabitto/Desktop/th17/mm_sub/mm_10k_53b_allpeaks_25sl_500s.bed'
f2name = '/Users/mgabitto/Desktop/th17/mm_sub/mm_5k_53b_allpeaks_25sl_500s.bed'
f3name = '/Users/mgabitto/Desktop/th17/mm_sub/mm_1k_53b_allpeaks_25sl_500s.bed'
f4name = '/Users/mgabitto/Desktop/th17/mm_sub/mm_500_53b_allpeaks_25sl_500s.bed'
bed53mm.append(bed_file_correlation(f1name, f2name, f3name, f4name, mouse_chromlens))

res53mm = np.zeros((len(bed53mm), 4))
for i in np.arange(len(bed53mm)):
    mmcor = []
    for k_ in bed53mm[i].keys():
        if not np.any(np.isnan(bed53mm[i][k_])):
            mmcor.append(bed53mm[i][k_])
    res53mm[i, :] = np.mean(np.array(mmcor), axis=0)[:, 0]
np.save('results/bed53_correl_exploremm.npy', res53mm)



bed75mm = []
f1name = '/Users/mgabitto/Desktop/th17/mm_sub/mm_10k_75b_allpeaks_100sl_250s.bed'
f2name = '/Users/mgabitto/Desktop/th17/mm_sub/mm_5k_75b_allpeaks_100sl_250s.bed'
f3name = '/Users/mgabitto/Desktop/th17/mm_sub/mm_1k_75b_allpeaks_100sl_250s.bed'
f4name = '/Users/mgabitto/Desktop/th17/mm_sub/mm_500_75b_allpeaks_100sl_250s.bed'
bed75mm.append(bed_file_correlation(f1name, f2name, f3name, f4name, mouse_chromlens))

f1name = '/Users/mgabitto/Desktop/th17/mm_sub/mm_10k_75b_allpeaks_100sl_500s.bed'
f2name = '/Users/mgabitto/Desktop/th17/mm_sub/mm_5k_75b_allpeaks_100sl_500s.bed'
f3name = '/Users/mgabitto/Desktop/th17/mm_sub/mm_1k_75b_allpeaks_100sl_500s.bed'
f4name = '/Users/mgabitto/Desktop/th17/mm_sub/mm_500_75b_allpeaks_100sl_500s.bed'
bed75mm.append(bed_file_correlation(f1name, f2name, f3name, f4name, mouse_chromlens))

f1name = '/Users/mgabitto/Desktop/th17/mm_sub/mm_10k_75b_allpeaks_50sl_250s.bed'
f2name = '/Users/mgabitto/Desktop/th17/mm_sub/mm_5k_75b_allpeaks_50sl_250s.bed'
f3name = '/Users/mgabitto/Desktop/th17/mm_sub/mm_1k_75b_allpeaks_50sl_250s.bed'
f4name = '/Users/mgabitto/Desktop/th17/mm_sub/mm_500_75b_allpeaks_50sl_250s.bed'
bed75mm.append(bed_file_correlation(f1name, f2name, f3name, f4name, mouse_chromlens))

f1name = '/Users/mgabitto/Desktop/th17/mm_sub/mm_10k_75b_allpeaks_50sl_500s.bed'
f2name = '/Users/mgabitto/Desktop/th17/mm_sub/mm_5k_75b_allpeaks_50sl_500s.bed'
f3name = '/Users/mgabitto/Desktop/th17/mm_sub/mm_1k_75b_allpeaks_50sl_500s.bed'
f4name = '/Users/mgabitto/Desktop/th17/mm_sub/mm_500_75b_allpeaks_50sl_500s.bed'
bed75mm.append(bed_file_correlation(f1name, f2name, f3name, f4name, mouse_chromlens))

f1name = '/Users/mgabitto/Desktop/th17/mm_sub/mm_10k_75b_allpeaks_25sl_250s.bed'
f2name = '/Users/mgabitto/Desktop/th17/mm_sub/mm_5k_75b_allpeaks_25sl_250s.bed'
f3name = '/Users/mgabitto/Desktop/th17/mm_sub/mm_1k_75b_allpeaks_25sl_250s.bed'
f4name = '/Users/mgabitto/Desktop/th17/mm_sub/mm_500_75b_allpeaks_25sl_250s.bed'
bed75mm.append(bed_file_correlation(f1name, f2name, f3name, f4name, mouse_chromlens))

f1name = '/Users/mgabitto/Desktop/th17/mm_sub/mm_10k_75b_allpeaks_25sl_500s.bed'
f2name = '/Users/mgabitto/Desktop/th17/mm_sub/mm_5k_75b_allpeaks_25sl_500s.bed'
f3name = '/Users/mgabitto/Desktop/th17/mm_sub/mm_1k_75b_allpeaks_25sl_500s.bed'
f4name = '/Users/mgabitto/Desktop/th17/mm_sub/mm_500_75b_allpeaks_25sl_500s.bed'
bed75mm.append(bed_file_correlation(f1name, f2name, f3name, f4name, mouse_chromlens))

res75mm = np.zeros((len(bed75mm), 4))
for i in np.arange(len(bed75mm)):
    mmcor = []
    for k_ in bed75mm[i].keys():
        if not np.any(np.isnan(bed75mm[i][k_])):
            mmcor.append(bed75mm[i][k_])
    res75mm[i, :] = np.mean(np.array(mmcor), axis=0)[:, 0]
np.save('results/bed75_correl_exploremm.npy', res75mm)








# ###############################################################
# ###############################################################
# ###############################################################
# ###############################################################
# ###############################################################
# ###############################################################
# PLOTS
# PLOTS
# PLOTS
# PLOTS


ll = []
ll.append('/Users/mgabitto/Dropbox (Simons Foundation)/2_PeakCalling/ATACseq_Data/10x_scTest/GM12878_mix_mm10/hg19/atac_v1_hg_10k_peaks.bed')
ll.append('/Users/mgabitto/Dropbox (Simons Foundation)/2_PeakCalling/ATACseq_Data/10x_scTest/GM12878_mix_mm10/hg19-5k/atac_v1_hg_5k_peaks.bed')
ll.append('/Users/mgabitto/Dropbox (Simons Foundation)/2_PeakCalling/ATACseq_Data/10x_scTest/GM12878_mix_mm10/hg19-1k/atac_v1_hg_1k_peaks.bed')
ll.append('/Users/mgabitto/Dropbox (Simons Foundation)/2_PeakCalling/ATACseq_Data/10x_scTest/GM12878_mix_mm10/hg19-500/atac_v1_hg_500_peaks.bed')
ll.append('/Users/mgabitto/Desktop/PD/hg19_10k_roundedtab.pd')
ll.append('/Users/mgabitto/Desktop/PD/hg19_5k_roundedtab.pd')
ll.append('/Users/mgabitto/Desktop/PD/hg19_1k_roundedtab.pd')
ll.append('/Users/mgabitto/Desktop/PD/hg19_500_roundedtab.pd')
ll.append('/Users/mgabitto/Desktop/MACS/10k/atac_v1_hg_10k_possorted_bam_peaks.bed')
ll.append('/Users/mgabitto/Desktop/MACS/5k/atac_v1_hg_5k_possorted_bam_peaks.bed')
ll.append('/Users/mgabitto/Desktop/MACS/1k/atac_v1_hg_1k_possorted_bam_peaks.bed')
ll.append('/Users/mgabitto/Desktop/MACS/500/atac_v1_hg_500_possorted_bam_peaks.bed')
ll.append('/Users/mgabitto/Desktop/th17/mm_sub/mm_10k_53b_allpeaks_100sl_500s.bed')
ll.append('/Users/mgabitto/Desktop/th17/mm_sub/mm_5k_53b_allpeaks_100sl_500s.bed')
ll.append('/Users/mgabitto/Desktop/th17/mm_sub/mm_1k_53b_allpeaks_100sl_500s.bed')
ll.append('/Users/mgabitto/Desktop/th17/mm_sub/mm_500_53b_allpeaks_100sl_500s.bed')

peaks = ll

chromlens = dh.human_lens()
chromlens.__delitem__('chr22')
chromlens.__delitem__('chrM')
chromlens.__delitem__('chrX')
chromlens.__delitem__('chrY')
length = []
frac = []
for p_ in peaks:
    bd = bt(p_)
    length.append(len(bd))

    cc = {}
    for key in chromlens:
        cc[key] = 0
    for interval in bd:
        if interval[0] in chromlens.keys():
            cc[interval[0]] += np.abs(np.int(interval[2]) - np.int(interval[1]))

    temp = []
    for key in chromlens:
        temp.append(np.float(cc[key]) / chromlens[key])
    frac.append(np.mean(temp))

pd=np.load("/Users/mgabitto/Desktop/Projects/PeakCalling/ChromA/results/count_chr1_21_PD.npy")
m=np.load("/Users/mgabitto/Desktop/Projects/PeakCalling/ChromA/results/count_chr1_21_M.npy")
chroma=np.load("/Users/mgabitto/Desktop/Projects/PeakCalling/ChromA/results/count_chr1_53_100_500.npy")

fig401 = plt.figure()
ax = fig401.subplots(3, 1)
cc = [pd[0], m[0], chroma[0],pd[1], m[1], chroma[1],pd[2], m[2], chroma[2], pd[3], m[3], chroma[3]]
ax[0].bar(['PB', 'M', 'CHROMA','PB1', 'M1', 'CHROMA1','PB2', 'M2', 'CHROMA2','PB3', 'M3', 'CHROMA3'], cc, width=0.2)
plt.xlabel('Algorithm')
plt.ylabel('Chipseq Recalled')

cc = [length[4], length[8], length[12], length[5], length[9], length[13], length[6], length[10], length[14], length[7], length[11], length[15],]
ax[1].bar(['PB', 'M', 'CHROMA','PB1', 'M1', 'CHROMA1','PB2', 'M2', 'CHROMA2','PB3', 'M3', 'CHROMA3'], cc, width=0.2)
plt.xlabel('Algorithm')
plt.ylabel('Number of peaks')

cc = [frac[4], frac[8], frac[12], frac[5], frac[9], frac[13], frac[6], frac[10], frac[14], frac[7], frac[11], frac[15],]
ax[2].bar(['PB', 'M', 'CHROMA','PB1', 'M1', 'CHROMA1','PB2', 'M2', 'CHROMA2','PB3', 'M3', 'CHROMA3'], cc, width=0.2)
plt.xlabel('Algorithm')
plt.ylabel('Accessible Genome Fraction')

pp = PdfPages('sc_cs_calls.pdf')
pp.savefig(fig401)
pp.close()






fig402 = plt.figure()
plt.plot(M_mm, '-+')
plt.plot(PD_mm, '-+')
plt.plot(bed53_mmcorrel[0, :], '-+')
plt.xlabel('Algorithm')
plt.ylabel('Chipseq Recalled')
plt.axis([-0.15, 3.15, 0.485, 1.15])
pp = PdfPages('sc_correlationMouse.pdf')
pp.savefig(fig402)
pp.close()


fig402 = plt.figure()
plt.plot(M_hg, '-+')
plt.plot(PD_hg, '-+')
plt.plot(bed53_hgcorrel[0, :], '-+')
plt.xlabel('Algorithm')
plt.ylabel('Chipseq Recalled')
plt.axis([-0.15, 3.15, 0.485, 1.15])
pp = PdfPages('sc_correlationHuman.pdf')
pp.savefig(fig402)
pp.close()


M_mm = []
for i in bedM_mmcorrel.keys():
    M_mm.append(bedM_mmcorrel[i])
M_mm = np.array(M_mm).mean(axis=0)
PD_mm = []
for i in bedPD_mmcorrel.keys():
    PD_mm.append(bedPD_mmcorrel[i])
PD_mm = np.array(PD_mm).mean(axis=0)


M_hg = []
for i in bedM_hgcorrel.keys():
    M_hg.append(bedM_hgcorrel[i])
M_hg = np.array(M_hg).mean(axis=0)
PD_hg = []
for i in bedPD_hgcorrel.keys():
    PD_hg.append(bedPD_hgcorrel[i])
PD_hg = np.array(PD_hg).mean(axis=0)


