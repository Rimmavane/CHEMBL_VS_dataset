import pandas as pd
from plotly.subplots import make_subplots
import plotly.graph_objects as go


def similarity_plotting(csv_path, title='Similarity heatmap', output_name='similarity_heatmap', threshold=0.4, size=None, font_size=dict()):
    main_table = pd.read_csv(csv_path, index_col=0)
    main_table = main_table.loc[(main_table > threshold).any(1)]
    print(main_table)
    fig = go.Figure(data=go.Heatmap(
        z=main_table.values.tolist(),
        x=main_table.index.tolist(),
        y=main_table.columns.tolist()))

    fig.update_layout(title=go.layout.Title(text=title, xref="paper", x=0.5))
    if size is not None:
        fig.update_layout(autosize=False, width=size[0], height=size[1])
    fig.update_layout(
        title=go.layout.Title(text=title, xref="paper", x=0.5, font={'size': int(font_size['size'] * 1.3)}),
        font=font_size)
    fig.show()
    fig.write_html(output_name + '.html')


def activities_plotting(csv_path, group, title, output_name, target_col='Name', targets_per_plot=40,
                        names=None, colors=None, group_legend=None, size=None, share_y=None, font_size=dict(),
                        y_title=False):
    """
    Make a relative bar plot for selected colmuns.
    """
    # INITILIZE
    main_table = pd.read_csv(csv_path, index_col=0)
    main_table = main_table[main_table['In_PDBBIND'] == 0][
        ['ChEMBL ID', 'Name', 'Kd_active', 'Kd_inactive', 'Ki_active',
         'Ki_inactive', 'IC50_active', 'IC50_inactive',
         'Active_compounds', 'Inactive_compounds']]
    data = main_table.to_dict(orient='list')
    num_of_subplots = int(len(data['ChEMBL ID']) / targets_per_plot) + (len(data['ChEMBL ID']) % targets_per_plot > 0)

    # PREPARE PLOT BASE AND DATA FOR PLOTTING
    # shared_yaxes - 'all' for shared Y axis for all subplots
    fig = make_subplots(rows=num_of_subplots, cols=1, row_heights=[1600 / num_of_subplots] * num_of_subplots,
                        shared_yaxes=share_y)

    targets_ids_names = [data[target_col][i:i + targets_per_plot] for i in
                         range(0, len(data[target_col]), targets_per_plot)]
    targets_ids = [data['ChEMBL ID'][i:i + targets_per_plot] for i in range(0, len(data[target_col]), targets_per_plot)]
    for batch in range(len(targets_ids)):
        for name in range(len(targets_ids[batch])):
            splitted_name = targets_ids_names[batch][name][:15]
            splitted_id = targets_ids[batch][name].replace('CHEMBL', '')
            targets_ids[batch][name] = ".".join([splitted_name] + [str(splitted_id)])
    values = dict()

    # split data to batches
    for member in range(len(group)):
        splited_values = [data[group[member]][i:i + targets_per_plot] for i in
                          range(0, len(data[group[member]]), targets_per_plot)]
        values[group[member]] = splited_values

    if names is None:
        names = group

    # assert lengths are equal
    assert len(group) == len(colors)
    assert len(group) == len(names)

    # feed each subplot with portion of values
    leg = True
    for batch in range(len(targets_ids)):
        # make stack for each target in batch
        for val_index in range(len(group)):
            grouping_by_legend = names[val_index] if group_legend is None else group_legend[val_index]
            fig.append_trace(go.Bar(x=targets_ids[batch], y=values[group[val_index]][batch], name=names[val_index],
                                    marker_color=None if colors is None else colors[val_index],
                                    legendgroup=grouping_by_legend, showlegend=leg), row=batch + 1, col=1)
            if batch == 0:
                leg = False
    if y_title:
        for sub_plot in range(len(targets_ids)):
            fig.update_yaxes(title_text=y_title, col=sub_plot)
    fig.update_layout(barmode='relative', title=go.layout.Title(text=title, xref="paper", x=0.5))
    if size is not None:
        fig.update_layout(autosize=False, width=size[0], height=size[1])
    fig.update_layout(
        title=go.layout.Title(text=title, xref="paper", x=0.5, font={'size': int(font_size['size'] * 1.3)}),
        font=font_size,
        yaxis=dict(tickformat=".2%"))
    fig.write_html(output_name + '.html')


def blast_plotting(blast_csv, group, title, output_name, threshold=45, targets_per_plot=20, share_y=False, names=None,
                   colors=None, size=None, font_size=None, y_title='', group_legend=None, sort_by=None):
    if font_size is None:
        font_size = {'size': 22}
    blast_results = pd.read_csv(blast_csv, index_col=0)
    blast_results = blast_results.fillna(0)  # fill missing values with 0
    blast_results = blast_results
    blast_results = blast_results[group]

    to_drop = []
    for ix, row in blast_results.iterrows():
        if max(row) < threshold:
            to_drop.append(ix)
    blast_results = blast_results.drop(to_drop).divide(100)
    if sort_by is not None:
        blast_results = blast_results.sort_values(sort_by, axis=0, ascending=False)

    data = blast_results.to_dict(orient='list')
    num_of_subplots = int(len(blast_results) / targets_per_plot) + (len(blast_results) % targets_per_plot > 0)

    fig = make_subplots(rows=num_of_subplots, cols=1, row_heights=[1600 / num_of_subplots] * num_of_subplots,
                        shared_yaxes=share_y)

    targets_ids = blast_results.index.to_list()
    print(len(targets_ids))
    targets_ids = [targets_ids[i:i + targets_per_plot] for i in range(0, len(blast_results), targets_per_plot)]

    values = dict()

    # split data to batches
    for member in range(len(group)):
        splited_values = [data[group[member]][i:i + targets_per_plot] for i in
                          range(0, len(data[group[member]]), targets_per_plot)]
        values[group[member]] = splited_values

    if names is None:
        names = group

    for batch in range(len(targets_ids)):
        leg = False
        if batch == 0:
            leg = True
        # make stack for each target in batch
        for val_index in range(len(group)):
            grouping_by_legend = names[val_index] if group_legend is None else group_legend[val_index]
            fig.append_trace(go.Bar(x=targets_ids[batch], y=values[group[val_index]][batch], name=names[val_index],
                                    marker_color=None if colors is None else colors[val_index],
                                    legendgroup=grouping_by_legend, showlegend=leg), row=batch + 1, col=1)

    if y_title != '':
        for sub_plot in range(len(targets_ids)):
            fig.update_yaxes(title_text=y_title, col=sub_plot, tickformat="%")
    fig.update_layout(barmode='group', title=go.layout.Title(text=title, xref="paper", x=0.5))
    if size is not None:
        fig.update_layout(autosize=False, width=size[0], height=size[1])
    fig.update_layout(
        title=go.layout.Title(text=title, xref="paper", x=0.5, font={'size': int(font_size['size'] * 1.3)}),
        font=font_size)

    fig.write_html(output_name + '.html')



blast_plotting(blast_csv='chembl_blast_results.csv',
               group=['identity%_dekois', 'identity%_dude'],
               title='Biggest identity score between ChEMBL targets and DEKOIS and DUDE entries.',
               output_name='blast_chembl_dekois_dude_sorted_by_dekois',
               threshold=50,
               targets_per_plot=21,
               share_y=True,
               size=(2000, 3000),
               font_size=dict(size=22),
               names=['DEKOIS identity %', 'DUDE identity %'],
               colors=['red', 'blue'],
               y_title='Identity percent',
               sort_by='identity%_dekois')

blast_plotting(blast_csv='chembl_blast_results.csv',
               group=['identity%_dekois', 'identity%_dude'],
               title='Biggest identity score between ChEMBL targets and DEKOIS and DUDE entries.',
               output_name='blast_chembl_dekois_dude_sorted_by_dude',
               threshold=50,
               targets_per_plot=21,
               share_y=True,
               size=(2000, 3000),
               font_size=dict(size=22),
               names=['DEKOIS identity %', 'DUDE identity %'],
               colors=['red', 'blue'],
               y_title='Identity percent',
               sort_by='identity%_dude')

# activities_plotting(csv_path='actives_number_sampling/targets_after_fingerprint_similarity5_tc0.95.csv',
#                     group=['Kd_active', 'Ki_active', 'IC50_active'],
#                     title='ChEMBL active compounds with Kd, Ki or IC50 standard value for each target.',
#                     output_name='actives_50', names=['Kd', 'Ki', 'IC50'], colors=['orange', 'blue', 'forestgreen'],
#                     font_size=dict(size=22),
#                     targets_per_plot=21,
#                     share_y='all',
#                     size=(2000, 3500),
#                     y_title="Number of compounds")

# activities_plotting(csv_path='actives_number_sampling/targets_after_fingerprint_similarity5_tc0.95.csv',
#                     group=['Ki_inactive', 'IC50_inactive', 'Kd_inactive'],
#                     title='ChEMBL inactive compounds with Kd, Ki or IC50 standard value for each target.',
#                     output_name='inactives_50', names=['Kd', 'Ki', 'IC50'], colors=['orange', 'blue', 'forestgreen'],
#                     font_size=dict(size=22),
#                     targets_per_plot=21,
#                     share_y='all',
#                     size=(2000, 3500),
#                     y_title="Number of compounds")
#
# activities_plotting(csv_path='actives_number_sampling/targets_after_fingerprint_similarity5_tc0.95.csv',
#                     group=['Active_compounds', 'Inactive_compounds'],
#                     title='ChEMBL all active and inactive compounds for each target.',
#                     output_name='both_50', names=['Active', 'Inactive'], colors=['blue', 'grey'],
#                     font_size=dict(size=22),
#                     targets_per_plot=21,
#                     share_y='all',
#                     size=(2000, 3500),
#                     y_title="Number of compounds")
