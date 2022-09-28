import numpy as np
import pandas as pd
import math
import plotly.graph_objects as go
import plotly.express as px
from plotly.subplots import make_subplots


class Analysis():

    def __init__(self, eeg, canoms, panoms):
        self.eeg = eeg
        self.panoms = panoms
        self.canoms = canoms
        self.plot_data = self.set_plot_data()

    @staticmethod
    def filt_plotting_data(df, alpha):
        """Filters plotting data by anomaly strength.

        Args:
            df (DataFrame): Plotting data.
            alpha (_type_): Strength threshold.

        Returns:
            DataFrame: Subset of df satisfying that anomalies be greater or
            equal than alpha.
        """
        df["anoms"] = np.where(df["strength"] <= alpha, np.nan, df["anoms"])
        return (df)

    @staticmethod
    def alpha_update(dfs, fig_data, alpha):
        """Updates plotting data of all channels
        in some plotly's figure data according to
        anomaly strength. 

        Args:
            dfs (list): List containing plot data of all channels.
            fig_data (Dictionary): Data attribute of a plotly Figure object.
            alpha (_type_): Anomaly threshold.
        Returns:
            list: A list of values to serve as updated plotting traces.
        """
        n = len(fig_data)
        anom_dfs = [Analysis.filt_plotting_data(
            x, alpha)["anoms"].tolist() for x in dfs]
        eegs = [x.iloc[:, 1].tolist() for x in dfs]
        result = [None] * n
        result[::2] = eegs
        result[1::2] = anom_dfs
        return (result)

    def get_max_strength(self, rounded=False):
        """Gets maximum anomaly strength in the analysis.

        Args:
            rounded (bool, optional): Round up the maximum value? Defaults to False.

        Returns:
            float: Maximum anomaly strength in the analysis.
        """
        s = max(self.panoms["strength"] + self.canoms["mean.change"])
        s = math.ceil(s / 10) * 10 if rounded else s
        return s

    def tsubs(self, s, e):
        """Time-based subsetting function that applies to the EEG data,
        the collective anomalies data and the point anomalies data
        simultaneously.

        Args:
            s (float): Starting time
            e (float): Ending time

        Returns:
            tuple: EEG data, collective and point anomalies data, 
            all subsetted from time s to e.
        """
        eeg = self.eeg.tsubs(s, e)
        canoms = self.canoms[self.canoms["_Time"].between(s, e)]
        panoms = self.panoms[self.panoms["_Time"].between(s, e)]

        return (eeg, canoms, panoms)

    def set_for_plot(self, channel, joint=True):
        """Structures analysis data for the purposes of plotting.

        Args:
            channel (int): Channel column index
            joint (bool, optional): Should collective and point anomalies be stored
            a in a single column? Defaults to True. If false, two different columns
            are created for the two anomaly types.

        Returns:
            DataFrame: A data frame containing a channel's EEG data
            and all info pertaining to its anomalies.

        """
        # Remove unwanted channels
        eeg = self.eeg.data.iloc[:, [0, channel]]
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
        """Plots an analyzed channel.

        Args:
            channel (int): Channel's column index.
            knob (bool, optional): Include a slider to filter anomalies by 
            strength? Defaults to True.
            views (bool, optional): Include buttons that frame the channel view
            to specific time-intervals? Defaults to True.
            marker_size (int, optional): Size of anomaly markers. Defaults to 5.

        Returns:
            Figure (plotly): A plotly figure representing a channel an its anomalies.
        """

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
                    args=[{'y': [Analysis.filt_plotting_data(df, i)["anoms"].tolist()]}, [
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
        """Sets plotting data of all channels. This method is called
        on initialization so as to perform the heavy computation only
        once.

        Returns:
            list: List containing plotting data of each EEG channel.
        """
        dfs = []
        for i in range(1, self.eeg.nchans + 1):
            df = self.set_for_plot(i)
            dfs.append(df)
        return dfs

    def plot_analysis(self):
        """Plots all anomalous channels.

        Returns:
            Figure (plotly): Plotly figure representing all
            analyzed channels.
        """

        fig = make_subplots(rows=self.eeg.nchans, cols=1,
                            shared_xaxes=True,
                            shared_yaxes=True,
                            row_titles=self.eeg.chans)

        # Add traces and store subplots plot data
        for i in range(1, self.eeg.nchans + 1):
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
                args=[{"y":  Analysis.alpha_update(
                    self.plot_data, fig.data, i)}],
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

        vis_list = [None] * self.eeg.nchans * 2
        vis_list[::2] = [True] * self.eeg.nchans
        vis_list[1::2] = [False] * self.eeg.nchans
        raw_vis_list = [True] * self.eeg.nchans * 2
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
