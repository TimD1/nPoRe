import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib.lines as lines

plt.rcParams.update({'font.size': 22})
prefix = "results/data/par-happy"
boxx = 0.2
boxw = 0.5

# beds = ['np_0', 'np_1', 'np_all', 'all']
# corners = [0.99, 0.97, 0.97, 0.99]
# zoomboth = [True, False, False, False]
# corner2 = 0.85
# evalbeds = ['eval41']
# callvcfs = ['g5', 'g5-hap', 'g5-ra-hap']
# colors = ['purple', 'red', 'green']
# truthvcfs = ['truth41', 'truth41std']
# markers = ['+', '.']

beds = ['all']
corners = [0.99]
zoomboth = [False]
corner2 = 0.85
evalbeds = ['eval41']
callvcfs = ['g5-orig1', 'g5-orig2']
colors = ['red', 'purple']
truthvcfs = ['truth41']
markers = ['.']
labels = ["clair3-pileup", "clair3-full"]


# new figure for each (bed,evalbed) region
for b, bed in enumerate(beds):
    for e, evalbed in enumerate(evalbeds):
        print(f"> plotting {bed}-{evalbed}")

        # create figure
        corner = corners[b]
        fig, axs = plt.subplots(1,2, figsize=(15,7))

        # new plot line for each (call,truth) set
        for c, callvcf in enumerate(callvcfs):
            for t, truthvcf in enumerate(truthvcfs):

                # plot SNP and INDELs
                files = [f"{prefix}/{callvcf}-{bed}-{truthvcf}-{evalbed}"
                                 ".roc.Locations.SNP.PASS.csv", 
                         f"{prefix}/{callvcf}-{bed}-{truthvcf}-{evalbed}"
                                 ".roc.Locations.INDEL.PASS.csv"]
                for f, filename in enumerate(files):
                    first = True
                    for line in open(filename, 'r'):
                        if first:
                            first = False
                            continue
                        fields = line.split(',')
                        try: # plot each data point
                            recall = float(fields[7])
                            prec = float(fields[8])
                            axs[f].plot(recall, prec, 
                                    color=colors[c], 
                                    marker=markers[t], 
                                    markerfacecolor=colors[c],
                                    linestyle='None')

                            # plot zoom box
                            if not f and recall > corners[b] and prec > corners[b]:
                                recall2 = boxx+boxw*(recall-corners[b])/(1-corners[b])
                                prec2 = boxx+boxw*(prec-corners[b])/(1-corners[b])
                                axs[0].plot(recall2, prec2, 
                                        color=colors[c], 
                                        marker=markers[t], 
                                        markerfacecolor=colors[c],
                                        linestyle='None')
                            elif f and zoomboth[b] and recall > corner2 and prec > corner2:
                                recall2 = boxx+boxw*(recall-corner2)/(1-corner2)
                                prec2 = boxx+boxw*(prec-corner2)/(1-corner2)
                                axs[1].plot(recall2, prec2, 
                                        color=colors[c], 
                                        marker=markers[t], 
                                        markerfacecolor=colors[c],
                                        linestyle='None')
                        except ValueError:
                            pass
                    axs[f].set_xlabel('Recall')
                    axs[f].set_ylabel('Precision')

        # add zoom squares
        axs[0].add_patch(patches.Rectangle((corner,corner),
                .999-corner,.999-corner,fill=False, linewidth=2))
        axs[0].plot([boxx, corner], [boxx+boxw, 1], color='k', linestyle=':')
        axs[0].plot([boxx+boxw, 1], [boxx, corner], color='k', linestyle=':')
        axs[0].add_patch(patches.Rectangle((boxx,boxx),
                boxw,boxw,fill=False, linewidth=2))
        if zoomboth[b]:
            axs[1].add_patch(patches.Rectangle((corner2,corner2),
                    .999-corner2,.999-corner2,fill=False, linewidth=2))
            axs[1].plot([boxx, corner2], [boxx+boxw, 1], color='k', linestyle=':')
            axs[1].plot([boxx+boxw, 1], [boxx, corner2], color='k', linestyle=':')
            axs[1].add_patch(patches.Rectangle((boxx,boxx),
                    boxw,boxw,fill=False, linewidth=2))

        # # get ref indel/snps
        # ref_summary_file = open(f"{prefix}/{callvcfs[0]}-{bed}-truth41-eval41.summary.csv", "r")
        # ref_summary = ref_summary_file.readlines()
        # ref_indels = ref_summary[2].split(',')[2]
        # ref_snps = ref_summary[4].split(',')[2]

        # # get std ref indel/snps
        # std_summary_file = open(f"{prefix}/{callvcfs[0]}-{bed}-truth41std-eval41.summary.csv", "r")
        # std_summary = std_summary_file.readlines()
        # std_indels = std_summary[2].split(',')[2]
        # std_snps = std_summary[4].split(',')[2]

        # # calculate title percentages
        # axs[0].set_title(f'SNPs:  $+${100*int(ref_snps)/162604:.2f}%'
        #                  f'  $\\bullet${100*int(std_snps)/159963:.2f}%')
        # axs[1].set_title(f'INDELs:  $+${100*int(ref_indels)/28841:.2f}%'
        #                  f'  $\\bullet${100*int(std_indels)/32200:.2f}%')

        # snp axes
        axs[0].set_xlim(0,1)
        axs[0].set_ylim(0,1)
        axs[0].set_xticks([0, 0.2, 0.4, 0.6, 0.8, 1])
        axs[0].set_xticklabels(['0%', '20%', '40%', '60%', '80%', '100%'])
        axs[0].set_yticks([0, 0.2, 0.4, 0.6, 0.8, 1])
        axs[0].set_yticklabels(['0%', '20%', '40%', '60%', '80%', '100%'])

        # indel axes
        axs[1].set_xlim(0,1)
        axs[1].set_ylim(0,1)
        axs[1].set_xticks([0, 0.2, 0.4, 0.6, 0.8, 1])
        axs[1].set_xticklabels(['0%', '20%', '40%', '60%', '80%', '100%'])
        axs[1].set_yticks([0, 0.2, 0.4, 0.6, 0.8, 1])
        axs[1].set_yticklabels(['0%', '20%', '40%', '60%', '80%', '100%'])

        # zoom box labels
        axs[0].text(boxx-0.01,boxx, f'{int(corners[b]*100)}%', ha='right', va='bottom')
        axs[0].text(boxx,boxx-0.01, f'{int(corners[b]*100)}%', ha='center', va='top')
        axs[0].text(boxx-0.01,boxx+boxw, f'100%', ha='right', va='center')
        axs[0].text(boxx+boxw,boxx-0.01, f'100%', ha='center', va='top')
        if zoomboth[b]:
            axs[1].text(boxx-0.01,boxx, f'{int(corner2*100)}%', ha='right', va='bottom')
            axs[1].text(boxx,boxx-0.01, f'{int(corner2*100)}%', ha='center', va='top')
            axs[1].text(boxx-0.01,boxx+boxw, f'100%', ha='right', va='center')
            axs[1].text(boxx+boxw,boxx-0.01, f'100%', ha='center', va='top')

        legend_elements = [ \
                patches.Patch(facecolor=c, edgecolor=c, label=l) \
                        for c, l in zip(colors, labels)]
        axs[1].legend(handles = legend_elements)

        # save
        plt.tight_layout()
        plt.savefig(f'results/img/happy_acc-{bed}-{evalbed}.png')
        plt.close()
