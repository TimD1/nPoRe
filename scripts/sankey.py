from collections import defaultdict

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import rc
rc('font', **{'family': 'DejaVu Serif', 'serif': ['Computer Modern']})
import numpy as np
import pandas as pd


def sankey(lefts, rights, colors, leftWeights=None, rightWeights=None,
           leftLabels=None, rightLabels=None, rightColors=None, gaps=None,
           fontsize=14, figureName=None, bottoms = None):
    '''
    Make Sankey Diagram showing flow from left-->right
    Inputs:
        left = NumPy array of object labels on the left of the diagram
        right = NumPy array of corresponding labels on the right of the diagram
            len(right) == len(left)
        leftWeight = NumPy array of weights for each strip starting from the
            left of the diagram, if not specified 1 is assigned
        rightWeight = NumPy array of weights for each strip starting from the
            right of the diagram, if not specified the corresponding leftWeight
            is assigned
        colorDict = Dictionary of colors to use for each label
            {'label':'color'}
        leftLabels = order of the left labels in the diagram
        rightLabels = order of the right labels in the diagram
        rightColor = If true, each strip in the diagram will be be colored
                    according to its left label
    Ouput:
        None
    '''

    # set default labels and weights
    if rightColors is None:
        rightColors = [False] * len(rights)
    assert(len(lefts) == len(rights) == len(rightWeights) == len(gaps) == \
            len(leftWeights) == len(rightColors) == len(bottoms)-1)

    plt.figure()
    plt.rc('text', usetex=False)
    plt.rc('font', family='serif')

    for idx, label in enumerate(bottoms):
        plt.text(
            idx, 0, label,
            {'ha': 'center', 'va': 'top'},
            fontsize=fontsize+2, fontweight='bold'
        )

    n = len(lefts)
    for idx, (left, right, leftWeight, rightWeight, rightColor) in \
            enumerate(zip(lefts, rights, leftWeights, rightWeights, rightColors)):

        # set weights
        if not len(leftWeight):
            leftWeight = np.ones(len(left))
        if not len(rightWeight):
            rightWeight = leftWeight

        # init df
        dataFrame = pd.DataFrame({'left': left, 'right': right, 
            'leftWeight': leftWeight, 'rightWeight': rightWeight}, 
            index=range(len(left)))

        print(idx)
        print(dataFrame)
        print(' ')

        # set labels
        leftLabels = pd.Series(dataFrame.left.unique()).unique()
        rightLabels = pd.Series(dataFrame.right.unique()).unique()

        # Determine widths of individual strips
        ns_l = defaultdict()
        ns_r = defaultdict()
        for leftLabel in leftLabels:
            leftDict = {}
            rightDict = {}
            for rightLabel in rightLabels:
                leftDict[rightLabel] = dataFrame[(dataFrame.left == leftLabel) & \
                        (dataFrame.right == rightLabel)].leftWeight.sum()
                rightDict[rightLabel] = dataFrame[(dataFrame.left == leftLabel) & \
                        (dataFrame.right == rightLabel)].rightWeight.sum()
            ns_l[leftLabel] = leftDict
            ns_r[leftLabel] = rightDict

        # Determine positions of left label patches and total widths
        leftWidths = defaultdict()
        for i, leftLabel in enumerate(leftLabels):
            myD = {}
            myD['left'] = dataFrame[dataFrame.left == leftLabel].leftWeight.sum()
            if i == 0:
                myD['bottom'] = 0.02 * dataFrame.leftWeight.sum()
                myD['top'] = myD['bottom'] + myD['left']
            else:
                myD['bottom'] = leftWidths[leftLabels[i - 1]]['top'] + \
                        0.02 * dataFrame.leftWeight.sum()
                myD['top'] = myD['bottom'] + myD['left']
                topEdge = myD['top']
            leftWidths[leftLabel] = myD

        # Determine positions of right label patches and total widths
        rightWidths = defaultdict()
        for i, rightLabel in enumerate(rightLabels):
            myD = {}
            myD['right'] = dataFrame[dataFrame.right == rightLabel].rightWeight.sum()
            if i == 0:
                myD['bottom'] = 0.02 * dataFrame.rightWeight.sum()
                myD['top'] = myD['bottom'] + myD['right']
                if myD['right']:
                    myD['top'] = myD['bottom'] + myD['right']
                else:
                    myD['top'] = 0
            else:
                myD['bottom'] = rightWidths[rightLabels[i - 1]]['top'] + \
                        0.02 * dataFrame.rightWeight.sum()
                myD['top'] = myD['bottom'] + myD['right']
                topEdge = myD['top']
            rightWidths[rightLabel] = myD

        # Total vertical extent of diagram
        bar_width = 0.01

        # Draw vertical bars on left and right of each  label's section & print label
        for leftLabel in leftLabels:
            if leftWidths[leftLabel]['left']: # skip if weight is zero
                if gaps[idx]:
                    plt.fill_between(
                        [idx-3*bar_width, idx-bar_width],
                        2 * [leftWidths[leftLabel]['bottom']],
                        2 * [leftWidths[leftLabel]['bottom'] + leftWidths[leftLabel]['left']],
                        color=colors[leftLabel],
                        alpha=0.99,
                        linewidth=0,
                    )
                    plt.fill_between(
                        [idx+bar_width, idx+3*bar_width],
                        2 * [leftWidths[leftLabel]['bottom']],
                        2 * [leftWidths[leftLabel]['bottom'] + leftWidths[leftLabel]['left']],
                        color=colors[leftLabel],
                        alpha=0.99,
                        linewidth=0,
                    )
                else:
                    plt.fill_between(
                        [idx-bar_width, idx+bar_width],
                        2 * [leftWidths[leftLabel]['bottom']],
                        2 * [leftWidths[leftLabel]['bottom'] + leftWidths[leftLabel]['left']],
                        color=colors[leftLabel],
                        alpha=0.99,
                        linewidth=0,
                    )
                if idx == 0:
                    ha = 'left'
                else:
                    ha = 'center'
                plt.text(
                    idx,
                    leftWidths[leftLabel]['bottom'] + 0.5 * leftWidths[leftLabel]['left'],
                    leftLabel,
                    {'ha': ha, 'va': 'center'},
                    fontsize=fontsize
                )
        if idx == n-1:
            for rightLabel in rightLabels:
                if rightWidths[rightLabel]['right']: # skip if weight is zero
                    plt.fill_between(
                        [idx+1-bar_width, idx+1+bar_width], 
                        2 * [rightWidths[rightLabel]['bottom']],
                        2 * [rightWidths[rightLabel]['bottom'] + rightWidths[rightLabel]['right']],
                        color=colors[rightLabel],
                        alpha=0.99,
                        linewidth=0,
                    )
                    plt.text(
                        idx+1,
                        rightWidths[rightLabel]['bottom'] + 0.5 * rightWidths[rightLabel]['right'],
                        rightLabel,
                        {'ha': 'right', 'va': 'center'},
                        fontsize=fontsize
                    )

        # Plot strips
        for leftLabel in leftLabels:
            for rightLabel in rightLabels:
                if leftWidths[leftLabel]['left'] and \
                        rightWidths[rightLabel]['right']: # skip if weight is zero
                    labelColor = leftLabel
                    if rightColor:
                        labelColor = rightLabel
                    if len(dataFrame[(dataFrame.left == leftLabel) & \
                            (dataFrame.right == rightLabel)]) > 0:
                        # Create array of y values for each strip, half at left value,
                        # half at right, convolve
                        ys_d = np.array(50 * [leftWidths[leftLabel]['bottom']] + \
                                50 * [rightWidths[rightLabel]['bottom']])
                        ys_d = np.convolve(ys_d, 0.05 * np.ones(20), mode='valid')
                        ys_d = np.convolve(ys_d, 0.05 * np.ones(20), mode='valid')
                        ys_u = np.array(50 * [leftWidths[leftLabel]['bottom'] + \
                                ns_l[leftLabel][rightLabel]] + \
                                50 * [rightWidths[rightLabel]['bottom'] + \
                                ns_r[leftLabel][rightLabel]])
                        ys_u = np.convolve(ys_u, 0.05 * np.ones(20), mode='valid')
                        ys_u = np.convolve(ys_u, 0.05 * np.ones(20), mode='valid')

                        # Update bottom edges at each label
                        leftWidths[leftLabel]['bottom'] += ns_l[leftLabel][rightLabel]
                        rightWidths[rightLabel]['bottom'] += ns_r[leftLabel][rightLabel]
                        start = idx+bar_width
                        end = idx+1-bar_width
                        if gaps[idx]:
                            start = idx + 3*bar_width
                        if idx < n-1 and gaps[idx+1]:
                            end = idx+1 - 3*bar_width

                        plt.fill_between(
                            np.linspace(start, end, \
                                    len(ys_d)), ys_d, ys_u, alpha=0.75,
                            facecolor='none', linewidth=0,
                            color=colors[labelColor]
                        )
    plt.gca().axis('off')
    plt.gcf().set_size_inches(3*n, 4.5)
    if figureName != None:
        plt.savefig("{}.png".format(figureName), bbox_inches='tight', dpi=300)
        plt.close()
