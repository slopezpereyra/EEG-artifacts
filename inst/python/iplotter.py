# Functions in this file are designed to deal with EEG data as formatted
# by the artifactor R package specifically.

import numpy as np
import pandas as pd
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from plotly_resampler import FigureResampler, register_plotly_resampler


def get_chan_data(df, chan):
    """Uses the linear function y = 3x - 2 to retrieve
    channel columns pertaining a specific channel x


    Args:
        df (DataFrame): An analysis dataframe as formatted by
        artifactor's set.iplot.data method.
        chan (int): Channel whose data will be retrieved.

    Returns:
        DataFrame: Subset of df containing only the columns 
        pertaining to channel chan.
    """
    s = chan * 3 - 2
    data = df.iloc[:, [0, s, s + 1, s + 2]]
    return data


def count_channels(df):
    """Counts how many channels are in the analysis
    data frame df.

    Args:
        df (DataFrame): An analysis dataframe as formatted by
        artifactor's set.iplot.data method.

    Returns:
        int: Number of EEG channels in df.
    """
    return int(len(df.columns[1:]) / 3)


def get_chan_names(df):
    """Gets name of EEG channels in df.

    Args:
        df (DataFrame): An analysis dataframe as formatted by
        artifactor's set.iplot.data method.

    Returns:
        list: List of EEG channel names in df.
    """
    indexes = []
    for i in range(1, len(df.columns[1:]) // 3 + 1):
        indexes.append(i * 3 - 3)
    return (df.columns[1:][indexes])


def get_max_strength(df):
    """Returns maximum anomaly strength in df.

    Args:
        df (DataFrame): An analysis dataframe as formatted by
        artifactor's set.iplot.data method.

    Returns:
        float: Maximum anomaly strength in the analysis.
    """
    strength_cols = [col for col in df if col.startswith('strength')]
    max = np.nanmax(np.nanmax(df[strength_cols].values, axis=1))
    return max


def filter_channel_by_alpha(df, alpha):
    """Filters a single channel by anomaly strength.

    Args:
        df (DataFrame): An channel's analysis dataframe as returned by
        the get_chan_data function.
        alpha (float): Strength threshold.

    Returns:
        DataFrame: Subset of df satisfying that anomalies be greater or
        equal than alpha.
    """
    df.iloc[:, 2] = np.where(df.iloc[:, 3] <= alpha, np.nan, df.iloc[:, 2])
    return df.iloc[:, 2]


def filter_data_by_alpha(df, alpha):
    """Filters an analysis data frame by anomaly strength.

    Args:
        df (DataFrame): An analysis dataframe as formatted by
        artifactor's set.iplot.data method.
        alpha (float): Strength threshold.

    Returns:
        DataFrame: Subset of df satisfying that anomalies be greater or
        equal than alpha.
    """
    strength_cols = [col for col in df if col.startswith('strength')]
    anom_cols = [col for col in df if col.startswith("anoms")]
    df[anom_cols] = np.where(df[strength_cols] >= alpha, df[anom_cols], np.nan)
    return (df)


def update_fig_by_alpha(df, fig_data, alpha):
    """Given a plotly figure depicting an artifactor analysis,
    get the y-coordinates of each of its traces when filtering
    anomalies by anomaly strength.

    Args:
        df DataFrame: The analysis dataframe used for the 
        analysis plot. It is always a dataframe of the form
        specified by artifactor's set.iplot.data method.
        fig_data (dictionary): The analysis plot figure's data.
        alpha (float): Strength threshold

    Returns:
        list: A list of lists. Each nth sublist contains the y
        coordinate of the nth trace configuring the analysis figure
        with anomalies filtered by strength.
    """
    n = len(fig_data)

    nchans = count_channels(df)

    filtered = filter_data_by_alpha(df, alpha)
    anom_dfs = []
    eegs = []
    for chan in range(1, nchans + 1):
        chan_data = get_chan_data(filtered, chan)
        anom_dfs.append(chan_data.iloc[:, 2].tolist())
        eegs.append(chan_data.iloc[:, 1].tolist())
    result = [None] * n
    result[::2] = eegs
    result[1::2] = anom_dfs
    return (result)


def plot_channel(df, channel, s=-1, e=-1, knob=True, views=True, marker_size=5, resample=False):
    """Plots an analyzed channel.

    Args:
        channel (int): Channel's column index.
        s (int): Initial second of the plot.
        e (int): Final second of the plot.
        knob (bool, optional): Include a slider to filter anomalies by
        strength? Defaults to True.
        views (bool, optional): Include buttons that frame the channel view
        to specific time-intervals? Defaults to True.
        marker_size (int, optional): Size of anomaly markers. Defaults to 5.
        resample(bool, optional): Resample the EEG for faster processing?

    Returns:
        Figure (plotly): A plotly figure representing a channel an its anomalies.
    """

    # Set plotting data
    df = get_chan_data(df, channel)
    if (s == -1):
        s, e = df["Time"].iloc[0], df["Time"].iloc[-1]
    df = df[df["Time"].between(s, e)]
    df["Time"] = pd.to_datetime(df["Time"] * 1000000000)
    # xy values for EEG and anoms.
    x0, y0 = df["Time"], df.iloc[:, 1]
    x1, y1 = df["Time"], df.iloc[:, 2]

    traces = [go.Scattergl(x=x0, y=y0, line=dict(width=1, color="black")),
              go.Scattergl(x=x1, y=y1, mode='markers', line=dict(color="red", width=0.8),
                           marker=dict(size=marker_size), visible=True, name="anoms")]
    if resample:
        fig = FigureResampler(go.Figure(data=traces))
    else:
        fig = go.Figure(data=traces)
    if knob:
        steps = []
        for i in range(0, round(np.nanmax(df.iloc[:, 3])) + 10, 10):
            step = dict(
                method='restyle',
                args=[{'y': [filter_channel_by_alpha(df, i).tolist()]}, [
                    1, 2]],
                label=i
            )
            steps.append(step)

        sliders = [dict(
            active=0,
            currentvalue={"prefix": "Mean change >= "},
            pad={"t": 200},
            steps=steps
        )]

        fig.update_layout(
            sliders=sliders
        )

    if views:
        fig.update_layout(
            updatemenus=[
                dict(
                    buttons=list([
                        dict(
                            args=["visible", True, [1, 2]],
                            label="Artifact view",
                            method="restyle"
                        ),
                        dict(
                            args=["visible", False, [1, 2]],
                            label="Raw view",
                            method="restyle"
                        ),
                    ]),
                    direction="down",
                    showactive=True,
                    xanchor="left",
                    yanchor="top"
                ),
            ]
        )
    fig.update_layout(showlegend=False)
    return (fig)


def plot_analysis(df, s=-1, e=-1, n_resample=2, marker_size=3, show=False, save=False, return_fig=False):
    """Plot all channels in a given analysis.

    Args:
        df (DataFrame): An analysis dataframe as formatted by
        artifactor's set.iplot.data method.
        s (int, optional): Initial second of the plot. Defaults to -1, which means
        the first second of the given analysis.
        e (int, optional): Last second of the plot. Defaults to -1,
        meaning the last second of the given analysis.
        n_resample (int, optional): Integer setting downsample decimation. If n,
        the only one every n values are plotted. Defaults to 2 and should be set
        to higher values for large analysis.
        marker_size (float, optional): Size of the red markers that signal
        anomalies in the plot. Defaults to 2.5.
        show (bool, optional): Show the plot immediately? Defaults to False.
        save (bool, optional): Save the plot in an .html file? Defaults to False.
        return_fig (bool, optional): Return the figure? Defaults to False.
        Since this function is designed to be used in the R package artifactor,
        it may be convenient to avoid lurking Python objects in the R environment.
        Therefore, the function should not return a plotly figure whenever used in
        the R package, and the conditional return clause is justified.

    Returns:
        _type_: _description_
    """

    df = df.iloc[::n_resample, :]
    nchans = count_channels(df)
    names = get_chan_names(df).tolist()
    fig = make_subplots(rows=nchans, cols=1,
                        shared_xaxes=True,
                        shared_yaxes=True,
                        row_titles=names,
                        vertical_spacing=0)

    # Add traces
    for i in range(1, nchans + 1):
        analysis = plot_channel(df, i, s, e, False, False,
                                marker_size=marker_size)
        fig.add_trace(analysis.data[0], row=i, col=1)
        fig.append_trace(analysis.data[1], row=i, col=1)

    # Add steps with updated subplot plot_data
    # based on alpha discrimination
    steps = []
    max = get_max_strength(df)
    # This may be wrong calculation of max
    for i in range(1, round(max) + 10, 10):
        step = dict(
            method='restyle',
            args=[{"y":  update_fig_by_alpha(
                df, fig.data, i)}],
            label=i
        )
        steps.append(step)

    sliders = [dict(
        active=0,
        currentvalue={"prefix": "Mean change >= "},
        pad={"t": 80},
        steps=steps,
        y=-0.1,
        yanchor="top"
    )]

    fig.update_layout(
        sliders=sliders
    )
    fig.update_layout(
        xaxis=dict(
            rangeselector=dict(
                buttons=list([
                    dict(count=1,
                         label="Minute",
                         step="minute",
                         stepmode="backward"),
                    dict(count=30,
                         label="Epoch",
                         step="second",
                         stepmode="backward"),
                    dict(count=1,
                         label="Hour",
                         step="hour",
                         stepmode="backward"),
                    dict(step="all")
                ])
            ),
            rangeslider=dict(
                visible=True
            ),
            type="date"
        )
    )
    for i in range(1, nchans):
        fig.update_xaxes(rangeslider={'visible': False}, row=i, col=1)
    fig.update_xaxes(rangeslider={'visible': True}, row=nchans, col=1)

    fig.update_layout(showlegend=False)
    if show:
        fig.show()
    if save:
        fig.write_html(
            "analysis_{}_to_{}_nresample={}.html".format(str(s), str(e), str(n_resample)))
    # This function is designed to be used in the R package artifactor. To avoid lurking
    # Python objects in the R environment, the function should not return a plotly figure
    # when called from R. Therefore, this conditional return logic is justified.
    if return_fig:
        return (fig)


def plot_eeg(df, s=-1, e=-1, joint=False, save=False, show=False):

    register_plotly_resampler(mode='auto', default_n_shown_samples=100000)

    if (s == -1 or e == -1):
        s, e = df["Time"].iloc[0], df["Time"].iloc[-1]

    nchans = len(df.columns) - 1
    names = df.columns[1:].tolist()
    df = df[df["Time"].between(s, e)]
    df["Time"] = pd.to_datetime(df["Time"] * 1000000000)
    x = df["Time"]

    fig = go.Figure() if joint else make_subplots(rows=nchans, cols=1,
                                                  shared_xaxes=True,
                                                  shared_yaxes=True,
                                                  row_titles=names,
                                                  vertical_spacing=0)

    if not joint:
        for i in range(1, nchans + 1):
            trace = go.Scattergl(x=x, y=df.iloc[:, i])
            fig.add_trace(trace, row=i, col=1)
    else:
        for i in range(1, nchans):
            trace = go.Scattergl(
                x=x, y=df.iloc[:, i], name=names[i - 1])
            fig.add_trace(trace)

    fig.update_layout(showlegend=joint)
    fig.update_annotations(font_size=10)

    # Add range slider
    fig.update_layout(
        xaxis=dict(
            rangeselector=dict(
                buttons=list([
                    dict(count=1,
                         label="Minute",
                         step="minute",
                         stepmode="backward"),
                    dict(count=30,
                         label="Epoch",
                         step="second",
                         stepmode="backward"),
                    dict(count=1,
                         label="Hour",
                         step="hour",
                         stepmode="backward"),
                    dict(step="all")
                ])
            ),
            rangeslider=dict(
                visible=True
            ),
            type="date"
        )
    )
    if not joint:
        for i in range(1, nchans):
            fig.update_xaxes(rangeslider={'visible': False}, row=i, col=1)
        fig.update_xaxes(
            rangeslider={'visible': True}, row=nchans, col=1)
    if save:
        joint_str = "joint" if joint else "sep"
        fig.write_html(
            "results/eeg_{}_to_{}_{}.html".format(str(s), str(e), joint_str))
    if show:
        fig.show()
    return fig
