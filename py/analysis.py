import numpy as np
import pandas as pd
import math
import plotly.graph_objects as go
import plotly.express as px
from plotly.subplots import make_subplots


class Analysis():

    def __init__(self, eeg, panoms, canoms):
        self.eeg = self._preprocess(eeg.copy())
        self.panoms = panoms.copy()
        self.canoms = canoms.copy()
        self.chans = self.eeg.columns.to_list()[1:]
        self.nchans = len(self.chans)
        self.plot_data = self.set_plot_data()

    def _preprocess(self, eeg):
        # The factor adjusts for datetime conversion from
        # data coming from the R artifactor package.
        eeg["Time"] = pd.to_datetime(eeg["Time"] * 1000000000)
        return (eeg)

    def get_max_strength(self, rounded=False):
        s = max(self.panoms["strength"] + self.canoms["mean.change"])
        s = math.ceil(s / 10) * 10 if rounded else s
        return s

    def set_for_plot(self, channel, joint=True):
        # Remove unwanted channels
        eeg = self.eeg.iloc[:, [0, channel]]
        canoms = self.canoms[self.canoms["variate"] == channel]
        panoms = self.panoms[self.panoms["variate"] == channel]

        # Insert empty anomaly and strength cols
        col_name = "anoms" if joint else "panoms"
        eeg.insert(2, col_name, np.nan)
        if not joint:
            eeg.insert(3, "canoms", np.nan)
        eeg.insert(len(eeg.columns), "strength", np.nan)

        # x location of each point anomaly
        points_x = panoms["location"].tolist()
        points_x = [x - 1 for x in points_x]
        # Fill anoms columns at anomalous x with the channels values
        eeg["anoms"].iloc[points_x] = eeg.iloc[:, 1].iloc[points_x]
        eeg["strength"].iloc[points_x] = panoms["strength"].tolist()

        col = "anoms" if joint else "canoms"
        for index, row in canoms.iterrows():
            s, e = int(row["start"]), int(row["end"])
            eeg[col][s:e] = eeg.iloc[:, 1][s:e]
            eeg["strength"][s:e] = row["mean.change"]

        return eeg

    def plot_channel(self, channel, knob=True, views=True, marker_size=5):

        # Set plotting data
        df = self.set_for_plot(channel)
        # xy values for EEG and anoms.
        x0, y0 = df["Time"], df.iloc[:, 1]
        x1, y1 = df["Time"], df.iloc[:, 2]

        traces = [go.Scatter(x=x0, y=y0, line=dict(color='black', width=1)),
                  go.Scatter(x=x1, y=y1, mode='lines+markers', line=dict(color="red", width=0.8),
                             marker=dict(size=marker_size), visible=True, name="anoms")]

        fig = go.Figure(data=traces)

        if knob:
            steps = []
            for i in range(0, round(np.nanmax(df["strength"])) + 2):
                step = dict(
                    method='restyle',
                    args=[{'y': [filt_plotting_data(df, i)["anoms"].tolist()]}, [
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
                        pad={"r": 10, "t": 10},
                        showactive=True,
                        x=0.1,
                        xanchor="left",
                        y=1.1,
                        yanchor="top"
                    ),
                ]
            )
        fig.update_layout(showlegend=False)
        return (fig)

    def set_plot_data(self):
        dfs = []
        for i in range(1, self.nchans + 1):
            df = self.set_for_plot(i)
            dfs.append(df)
        return dfs

    def plot_analysis(self):

        fig = make_subplots(rows=self.nchans, cols=1,
                            shared_xaxes=True,
                            shared_yaxes=True,
                            row_titles=self.chans)

        # Add traces and store subplots plot data
        for i in range(1, self.nchans + 1):
            analysis = self.plot_channel(i, False, False,
                                         marker_size=2.5)
            fig.add_trace(analysis.data[0], row=i, col=1)
            fig.append_trace(analysis.data[1], row=i, col=1)

        # Add steps with updated subplot plot_data
        # based on alpha discrimination
        steps = []
        for i in range(1, self.get_max_strength(True) + 1, 10):
            step = dict(
                method='restyle',
                args=[{"y":  alpha_update(self.plot_data, fig.data, i)}],
                label=i
            )
            steps.append(step)

        sliders = [dict(
            active=0,
            currentvalue={"prefix": "Mean change >= "},
            pad={"t": 0.1},
            steps=steps,
            y=-0.1,
            yanchor="top"
        )]

        fig.update_layout(
            sliders=sliders
        )

        vis_list = [None] * self.nchans * 2
        vis_list[::2] = [True] * self.nchans
        vis_list[1::2] = [False] * self.nchans
        raw_vis_list = [True] * self.nchans * 2
        fig.update_layout(
            updatemenus=[
                dict(
                    buttons=list([
                        dict(
                            args=["visible", raw_vis_list],
                            label="Artifact view",
                            method="restyle"
                        ),
                        dict(
                            args=["visible", vis_list],
                            label="Raw view",
                            method="restyle"
                        ),
                    ]),
                    direction="down",
                    pad={"r": 10, "t": 10},
                    showactive=True,
                    x=0.1,
                    xanchor="left",
                    y=1.1,
                    yanchor="top"
                ),
            ]
        )

        fig.update_layout(showlegend=False, height=1200)
        return (fig)


def filt_plotting_data(df, alpha):
    df["anoms"] = np.where(df["strength"] <= alpha, np.nan, df["anoms"])
    return (df)


def alpha_update(dfs, fig_data, alpha):
    n = len(fig_data)
    anom_dfs = [filt_plotting_data(x, alpha)["anoms"].tolist() for x in dfs]
    eegs = [x.iloc[:, 1].tolist() for x in dfs]
    result = [None] * n
    result[::2] = eegs
    result[1::2] = anom_dfs

    return (result)
