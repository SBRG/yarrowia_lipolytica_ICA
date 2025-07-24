from pymodulon.plotting import *
from pymodulon.util import _parse_sample, dima, explained_variance, mutual_info_distance
import matplotlib as mpl
mpl.rcParams['axes.linewidth'] = 1.5



def plot_dima_local(
    ica_data,
    sample1,
    sample2,
    threshold=5,
    fdr=0.1,
    label=True,
    adjust=True,
    table=False,
    color_inputs='blue',
    xlabel_to_plot=None,
    ylabel_to_plot=None,
    **kwargs,
):
    """
    Plots a DiMA plot between two projects or two sets of samples

    Parameters
    ----------
    ica_data: ~pymodulon.core.IcaData
        :class:`~pymodulon.core.IcaData` object
    sample1: list or str
        List of sample IDs or name of "project:condition" for x-axis
    sample2: list or str
        List of sample IDs or name of "project:condition" for y-axis
    threshold: float
        Minimum activity difference to determine DiMAs (default: 5)
    fdr: float
        False detection rate (default: 0.1)
    label: bool
        Label differentially activated iModulons (default: True)
    adjust: bool
        Automatically adjust labels (default: True)
    table: bool
        Return differential iModulon activity table (default: False)
    **kwargs:
        Additional keyword arguments passed to :func:`pymodulon.plotting.scatterplot`

    Returns
    -------
    ax: ~matplotlib.axes.Axes
        :class:`~matplotlib.axes.Axes` containing the scatterplot

    df_diff: ~pandas.DataFrame, optional
        Table reporting differentially activated iModulons

    """

    # use secret option to enable passing of clustered activity matrix
    try:
        A_to_use = kwargs.pop("alternate_A")
    except KeyError:
        A_to_use = ica_data.A

    # Override specific kwargs (their implementation is different
    # in this function)
    legend_cgw = kwargs.pop("legend", False)
    legend_kwargs_cgw = kwargs.pop("legend_kwargs", {})

    kwargs["legend"] = False
    kwargs["legend_kwargs"] = None

    # Get x and y coordinates
    sample1_list = _parse_sample(ica_data, sample1)
    sample2_list = _parse_sample(ica_data, sample2)
    if isinstance(sample1, str):
        xlabel = sample1
    else:
        xlabel = "\n".join(sample1)
    if isinstance(sample2, str):
        ylabel = sample2
    else:
        ylabel = "\n".join(sample2)

    a1 = A_to_use[sample1_list].mean(axis=1)
    a2 = A_to_use[sample2_list].mean(axis=1)

    df_diff = dima(
        ica_data,
        sample1_list,
        sample2_list,
        threshold=threshold,
        fdr=fdr,
        alternate_A=A_to_use,
    )

    groups = {}
    for i in A_to_use.index:
        if i not in df_diff.index:
            groups.update({i: "hidden"})
        else:
            groups.update({i: ""})

    scatter_kwargs = {
    's': 100,       # Size of scatter points
    'alpha': 1,   # Transparency
    'edgecolor': 'none',  # Edge color (if needed)
    }

    ax = scatterplot(
        a1,
        a2,
        groups=groups,
        line45=True,
        xlabel=xlabel_to_plot,
        ylabel=ylabel_to_plot,
        colors=color_inputs,  # <-- specify color here
        scatter_kwargs = scatter_kwargs,
        **kwargs,
    )

    # Get indices of hidden points (grey)
    hidden_indices = [i for i, g in groups.items() if g == "hidden"]
    
    # Update alpha for grey points
    for coll in ax.collections:
        if coll.get_offsets().shape[0] == len(hidden_indices):
            coll.set_alpha(0.2)  # Adjust the alpha value as needed

    if label:
        df_diff = pd.concat([df_diff, a1, a2], join="inner", axis=1)
        texts = []
        for k in df_diff.index:
            texts.append(ax.text(df_diff.loc[k, 0], df_diff.loc[k, 1], k, fontsize=9.8))
        if adjust:
            expand_args = {
                "expand_objects": (1.2, 1.4),
                "expand_points": (1.3, 1.3),
                "expand_text": (1.4, 1.4),
            }
            adjust_text(
                texts,
                ax=ax,
                arrowprops=dict(arrowstyle="->", color="grey", lw=1.5),
                only_move={"objects": "y"},
                **expand_args,
            )
    print(df_diff)
    # Add legend if requested
    if legend_cgw:
        ax.legend(**legend_kwargs_cgw)

    if table:
        return ax, df_diff

    else:
        return ax, df_diff
