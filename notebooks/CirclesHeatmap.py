'''
Created on 1/14/2018, 2018

@author: Ignacio Ibarra Del Rio

Description:
'''

class CirclesHeatmap():
    @staticmethod
    def plot(hm_values, hm_colors, xlabels, ylabels=None):
        if ylabels is None:
            ylabels = xlabels

        sns.set_style('ticks')

        coordinates = [[yi, xi] for yi, v in enumerate(hm_values) for xi in range(len(v))]
        values = [xi * 100 for yi, v in enumerate(hm_values) for xi in v]
        colors = [xi for yi, v in enumerate(hm_colors) for xi in v]
        # print coordinates
        # print values

        cm = plt.cm.get_cmap('RdBu')
        cm = truncate_colormap(cm, minval=-2.5, maxval=2.5)
        g = plt.scatter([p[0] for p in coordinates], [p[1] for p in coordinates],
                         c=colors, s=values, cmap=cm, vmin=-5, vmax=5,
                         edgecolors='black')
        plt.xticks([i for i in range(len(xlabels))], xlabels)
        plt.yticks([i for i in range(len(ylabels))], ylabels)

        plt.colorbar(g, label='Z-expression')

        labs = []
        g = []
        for size in [-5, -2.5, 0, 2.5, 5]:
            rgba = cm(size)

            print(size, rgba)
            g.append(plt.scatter([], [], s=200, marker='o', color=rgba[:3]))
            labs.append(str(size))
        plt.legend(g,
                   labs,
                   scatterpoints=1,
                   ncol=1,
                   fontsize=12,
                   bbox_to_anchor=(1.3, 1.05))
        plt.tight_layout()
        despine_all()

    @staticmethod
    def make_bubble_heatmap(order_frame, sizeDict, na_color='gray', title='title',
                            tickscolorbar=[-2, -1, 0, 1, 2], vmin=-2.5, vmax=2.5,
                            heatmap_grid=[2, 4, 0, 2, 2, 1],
                            circle_legend_grid=[2, 4, 0, 2, 2, 1],
                            colorbar_grid=[2, 5, 0, 3, 2, 1],
                            palette_id='RdBu_r', cbar_label='cbar_label', ncols=8,
                            marker=None, **kwargs):
        from matplotlib import gridspec
        sns.set_style('white')

        quantAmplifier = kwargs.get('quantAmplifier', 1) #factor for size of bubbles

        #function you received on the color gradient
        pallete = None
        if palette_id == 'RdBu_r':
            pallete = plt.cm.RdBu_r
        elif palette_id == 'Blues':
            pallete = plt.cm.Blues
        elif palette_id == 'Reds':
            pallete = plt.cm.Reds
        elif palette_id == 'Greens':
            pallete = plt.cm.Greens
        elif palette_id == 'Purples':
            pallete = plt.cm.Purples
        elif palette_id == 'YlGnBu':
            pallete = plt.cm.YlGnBu
        elif palette_id == 'PuOr':
            pallete = plt.cm.PuOr
        else:
            palette = palette_id
        scalarmap,colorList = CirclesHeatmap.get_specific_color_gradient(pallete,
                                                                         np.array(order_frame),
                                                                         vmin=vmin,
                                                                         vmax=vmax)

        fig = kwargs.get('fig', None)
        if fig is None:
            plt.clf()
            fig = plt.figure(figsize=(kwargs.get('w', 11), kwargs.get('h', 5)))


        nrows, ncols, rowi, coli, rowspan, colspan = heatmap_grid
        gs = gridspec.GridSpec(nrows, ncols)
        ax = plt.subplot(gs[rowi: rowi + rowspan, coli: coli + colspan])


        ylabelList = []
        i = 0

        # print order_frame
        # print sizeDict
        # to keep the order of dataframes in the same order as sizes, revert the dataframes
        order_frame = order_frame.reindex(index=order_frame.index[::-1])
        sizeDict = sizeDict.reindex(index=sizeDict.index[::-1])

        min_circle_size = kwargs.get('min_circle_size')
        max_circle_size = kwargs.get('max_circle_size')

        print('before plotting dimensions...')
        print(order_frame.shape)
        for ri, r in order_frame.iterrows():
            patList = list(r.values)

            sizeList = []
            if min_circle_size is not None and max_circle_size is not None:
                for si in list(sizeDict.loc[ri]):
                    if si < min_circle_size:
                        si = min_circle_size
                    scaled = (si - min_circle_size) / (max_circle_size - min_circle_size)
                    sizeList.append(scaled * quantAmplifier)
            else:
                sizeList = [(si ** kwargs.get('power', 1.0)) * (quantAmplifier) for si in list(sizeDict.loc[ri])]
            # print patList

            colorList = scalarmap.to_rgba(patList)
            colorList = [ci if not np.isnan(vi) else na_color
                         for ci, vi in zip(colorList, patList)]
            edgecolorList = list()
            x = list(range(len(patList)))
            y = [i] * len(patList)
            # print x, y, patList
            # print 'sizes'

            # define hataches, linewidths and alphas
            hatches = [None if np.isnan(pi) else None for pi in patList]
            linewidths = [0.1 if np.isnan(pi) else kwargs.get('sig_line_width', 2.0) for pi in patList]\
                if kwargs.get('line_widths', None) is None else list(kwargs.get('line_widths').loc[ri])
            # linewidths = [0.1 if np.isnan(pi) else 1.0 for pi in patList]
            alphas = [.5 if np.isnan(alpha) else 1.0 for alpha in patList]

            # print list(r.values)
            # print len(x), len(y), len(sizeList), len(colorList), len(hatches), len(linewidths), len(alphas)
            for xi, yi, si, ci, hatch_i,\
                    lw, alpha in zip(x, y, sizeList, colorList,
                              hatches, linewidths, alphas):
                # print xi, yi, si, ci
                # print 'power', kwargs.get('power', 1.0)
                # print(marker)
                ax.scatter(xi, yi, s=abs(si) * quantAmplifier,
                           marker=marker if marker is not None else ('o' if si < 0 else 'o'),
                           hatch=hatch_i, alpha=alpha,
                           edgecolor=kwargs.get('edgecolor', 'black'), color=ci, linewidth=lw)
                # print xi, yi, ri # list(order_frame.ix[ri])[si]
            ylabelList.append(ri)
            i += 1

        if kwargs.get('grid', True):
            plt.grid(True, linewidth=kwargs.get('grid_linewidth', 1.0))
        plt.title(kwargs.get('heatmap_title', 'title'))
        plt.xlabel(kwargs.get('xlab', 'xlab'), fontsize=14)
        plt.ylabel(kwargs.get('ylab', 'ylab'), fontsize=14)


        # print 'ylabels', ylabelList
        plt.yticks(list(range(len(ylabelList))))
        ax.set_yticklabels(ylabelList, fontsize=kwargs.get('yticks_fontsize', 11), color="black",
                           ha=kwargs.get('ha_ylabs', 'right'))
        # ax.set_xlim(-1,6)
        remove_top_n_right_ticks(ax)


        # print 'current columns'
        # print list(range(len(order_frame.columns)))
        # print order_frame.columns
        plt.xticks(list(range(len(order_frame.columns))), order_frame.columns,
                   fontsize=kwargs.get('xticks_fontsize', 11), rotation=kwargs.get('rotation_xlabs', 90), color="black",
                   ha=kwargs.get('ha_xlabs', 'center'))
        # ax.set_ylim(-2,len(ylabelList))

        lh, lt = get_legendHandle_for_second_sanity_check_plot(quantAmplifier=quantAmplifier * quantAmplifier,
                                                               marker='o' if marker is None else marker,
                                                               fmt=kwargs.get('fmt_legend', '%.1f'),
                                                               lw=0.5,
                                                               min_circle_size=min_circle_size,
                                                               max_circle_size=max_circle_size,
                                                               values=[tick ** kwargs.get('power', 1.0) for tick in kwargs.get('circle_legend_ticks')],
                                                               labels=kwargs.get('circle_legend_ticks'))
        if kwargs.get('show_circle_legend', True):
            l = plt.legend(lh, lt, bbox_to_anchor=kwargs.get('circle_legend_bbox', (1.8, 1)), scatterpoints=1,
                          title=kwargs.get('circles_legend_title', 'circles_legend_title'), ncol=1,
                          frameon=False)
            # Add the legend manually to the current Axes.
            ax = plt.gca().add_artist(l)

        # this is to add the circle (significant or not
        # print 'here...'
        if kwargs.get('show_sig_legend', False):
            # print kwargs.get('circle_legend_ticks')
            # print max(kwargs.get('circle_legend_ticks'))
            lh, lt = get_legendHandle_for_second_sanity_check_plot(quantAmplifier=quantAmplifier * quantAmplifier,
                                                                   marker='^' if marker is None else marker,
                                                                   fmt=kwargs.get('sig_fmt_legend', '%s'),
                                                                   lw=[kwargs.get('sig_line_width', 2.0), 0.0],
                                                                   labels=kwargs.get('sig_legend_ticks',
                                                                                      ['Yes', 'No']),
                                                                   values=kwargs.get('sig_legend_ticks',
                                                                                     ['Yes', 'No']),
                                                                   edgecolor=kwargs.get('edgecolor', 'black'),
                                                                   min_size_default=max(kwargs.get('circle_legend_ticks')) ** kwargs.get('power', 1.0))

            plt.legend(lh, lt, bbox_to_anchor=kwargs.get('sig_legend_bbox', (2.4, 1)),
                       title=kwargs.get('sig_legend_title', 'significant'), ncol=1, scatterpoints=1,
                       frameon=False)

        # ax1 = plt.subplot(gs[row_i + 2: row_i + 3, 4 + (plot_i * 8):6 + (plot_i * 8)])
        # ax1 = plt.subplot(gs[-2: -1, -(ncols / 2):-1])
        nrows, ncols, rowi, coli, rowspan, colspan = colorbar_grid
        ax1 = plt.subplot2grid([nrows, ncols], [rowi, coli], rowspan=rowspan, colspan=colspan)
        plt.axis('off')

        ax1.set_xticklabels([])
        ax1.set_yticklabels([])

        if kwargs.get('show_colorbar', True):
            # print 'ticks for colorbar:', tickscolorbar
            cbar = fig.colorbar(scalarmap, orientation="horizontal", format=kwargs.get('cbar_fmt_ticks', "%.1f"),
                                ticks = tickscolorbar)
            cbar.ax.tick_params(labelsize=kwargs.get('colorbar_ticks_labelsize', 12))
            cbar.set_label(cbar_label, fontsize=12)
        despine_all()

        # plt.tight_layout()

    @staticmethod
    def get_specific_color_gradient(colormap,inputList,**kwargs):
        vmin = kwargs.get('vmin','blaq')
        vmax = kwargs.get('vmax','blaq')
        cm = plt.get_cmap(colormap)
        if vmin=='blaq' or vmax=='blaq':
            if type(inputList)==list:
                cNorm = matplotlib.colors.Normalize(vmin=min(inputList), vmax=max(inputList))
            else:
                cNorm = matplotlib.colors.Normalize(vmin=inputList.min(), vmax=inputList.max())
        else:
            cNorm = matplotlib.colors.Normalize(vmin=vmin, vmax = vmax)
        scalarMap = matplotlib.cm.ScalarMappable(norm=cNorm, cmap=cm)
        scalarMap.set_array(inputList)
        colorList=scalarMap.to_rgba(inputList)
        return scalarMap,colorList


import networkx as nx
def repel_labels(ax, x, y, labels, k=0.01, arrowstyle="->", color='red', fontsize=8,
                 init_dy=0.0):
    G = nx.DiGraph()
    data_nodes = []
    init_pos = {}
    for xi, yi, label in zip(x, y, labels):
        data_str = 'data_{0}'.format(label)
        G.add_node(data_str)
        G.add_node(label)
        G.add_edge(label, data_str)
        data_nodes.append(data_str)
        init_pos[data_str] = (xi, yi)
        init_pos[label] = (xi, yi + init_dy)

    pos = nx.spring_layout(G, pos=init_pos, fixed=data_nodes, k=k)

    # undo spring_layout's rescaling
    pos_after = np.vstack([pos[d] for d in data_nodes])
    pos_before = np.vstack([init_pos[d] for d in data_nodes])
    scale, shift_x = np.polyfit(pos_after[:,0], pos_before[:,0], 1)
    scale, shift_y = np.polyfit(pos_after[:,1], pos_before[:,1], 1)
    shift = np.array([shift_x, shift_y])
    for key, val in pos.iteritems():
        pos[key] = (val*scale) + shift

    for label, data_str in G.edges():
        ax.annotate(label,
                    xy=pos[data_str], xycoords='data', fontsize=fontsize,
                    xytext=pos[label], textcoords='data',
                    arrowprops=dict(arrowstyle=arrowstyle,
                                    shrinkA=0, shrinkB=0,
                                    connectionstyle="arc3",
                                    color=color), )