import rpy2.robjects as ro
import pandas as pd
from py.analysis import Analysis
from rpy2.robjects import pandas2ri
from rpy2.robjects.conversion import localconverter
import plotly.graph_objects as go
import plotly.express as px
from plotly.subplots import make_subplots
import os


ro.r("""source('R/artifactor.R')""")


def load_eeg(data, signals, na_omit=False, exclude=[]):
    signals = pd.read_csv(signals)
    chan_names = signals["Label"].to_list()
    df = pd.read_csv(data, header=0,
                     names=["Time"] + chan_names)
    if na_omit:
        df = df.dropna()

    if len(exclude) > 0:
        df = df.drop(df.columns[exclude], axis=1)

    return df


class EEG():

    def __init__(self, data):
        self.data = data
        self.chans = self.data.columns.to_list()[1:]
        self.nchans = len(self.chans)

    def tsubs(self, s, e):
        sub = self.data[self.data["Time"].between(s, e)]
        return sub

    def analyze(self, s, e):
        segm = self.tsubs(s, e)
        # Convert to R data frame
        with localconverter(ro.default_converter + pandas2ri.converter):
            r_eeg = ro.conversion.py2rpy(segm)

        # Create r eeg object
        r_eeg = ro.r.new("eeg", data=r_eeg)
        # Create r analysis object
        analysis = ro.r.analyze(r_eeg, s, e, alpha=1)
        canoms = analysis.slots["canoms"]
        panoms = analysis.slots["panoms"]

        with localconverter(ro.default_converter + pandas2ri.converter):
            canoms = ro.conversion.rpy2py(canoms)
            panoms = ro.conversion.rpy2py(panoms)

        analysis = Analysis(segm, panoms, canoms)
        return (analysis)

    def draw(self, s, e, joint=False, width=1000, height=600):

        data = self.tsubs(s, e)
        data["Time"] = pd.to_datetime(data["Time"] * 1000000000)
        x = data["Time"]

        fig = go.Figure() if joint else make_subplots(rows=self.nchans, cols=1,
                                                      shared_xaxes=True,
                                                      shared_yaxes=True,
                                                      row_titles=self.chans)

        if not joint:
            for i in range(1, self.nchans + 1):
                trace = go.Scatter(x=x, y=data.iloc[:, i])
                fig.add_trace(trace, row=i, col=1)
        else:
            for i in range(1, self.nchans):
                trace = go.Scatter(
                    x=x, y=data.iloc[:, i], name=self.chans[i - 1])
                fig.add_trace(trace)

        fig.update_layout(showlegend=joint, width=width, height=height)
        fig.update_annotations(font_size=10)

        # Add range slider
        fig.update_layout(
            xaxis=dict(
                rangeselector=dict(
                    buttons=list([
                        dict(count=1,
                             label="Minute",
                             step="minute",
                             stepmode="todate"),
                        dict(count=60,
                             label="Hour",
                             step="minute",
                             stepmode="todate"),
                        dict(count=30,
                             label="Epoch",
                             step="second",
                             stepmode="todate"),
                        dict(step="all")
                    ])
                ),
                rangeslider=dict(
                    visible=False,

                ),
                type="date"
            )
        )
        return fig
