from urllib.request import DataHandler
import rpy2.robjects as ro
import pandas as pd
from analysis import Analysis
from rpy2.robjects import pandas2ri
from rpy2.robjects.conversion import localconverter
import plotly.graph_objects as go
import plotly.express as px
from plotly.subplots import make_subplots
import os

os.chdir("/home/santi/work/EEG-artifacts")
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
    """EEG Class
    """

    def __init__(self, data):
        """
        Args:
            data (Pandas DataFrame): DataFrame containing
            an EEG reading.
        """
        self.data = data
        self.chans = self.data.columns.to_list()[1:]
        self.nchans = len(self.chans)

    def tsubs(self, s, e):
        """Time-based subsetting function.

        Args:
            s (float): Starting time
            e (float): Ending time

        Returns:
            DataFrame: EEG data subsetted from second s
            to second e. 
        """
        sub = self.data[self.data["Time"].between(s, e)]
        return sub

    def _anoms_to_pandas(self, r_analysis):
        """Helper function that takes an R analysis object
        from the artifactor package and returns the @canoms
        and @panoms data frames converted to pandas format.

        Args:
            r_analysis (rObject): An R analysis object

        Returns:
            Tupple: Pair of pandas DataFrames of the form 
            (canoms, panoms) where canoms contains collective
            anomalies and panoms point anomalies from an
            EEG analysis.
        """
        with localconverter(ro.default_converter + pandas2ri.converter):
            canoms = ro.conversion.rpy2py(r_analysis.slots["canoms"])
            panoms = ro.conversion.rpy2py(r_analysis.slots["panoms"])

        return canoms, panoms

    def _rartifactor_analysis(self, s, e):
        """Run M-CAPA analysis over a subset of the
        EEG data by calling the R artifactor package

        Args:
            s (float): Starting time
            e (float): Ending time

        Returns:
            rObject: An rObject of class analysis (see artifactor package).
        """
        segm = self.tsubs(s, e)
        # Convert to R data frame
        with localconverter(ro.default_converter + pandas2ri.converter):
            r_eeg = ro.conversion.py2rpy(segm)

        # Create r eeg object
        r_eeg = ro.r.new("eeg", data=r_eeg)
        # Create r analysis object
        analysis = ro.r.analyze(r_eeg, s, e, alpha=1)
        return analysis

    def analyze(self, s, e):
        """Performs M-CAPA analysis, transform results from R
        data frames to pandas format and creates an instance of
        the Analysis class with them.

        Args:
            s (float): Starting time
            e (float): Ending time

        Returns:
            Analysis: An instance of class Analysis whose data is the
            subset of the EEG upon which M-CAPA analysis was performed,
            and whose canoms and panoms attributes are the collective
            and point anomalies found, respectively.
        """

        # R analysis object
        r_analysis = self._rartifactor_analysis(s, e)
        # canoms and anoms data in pandas form
        anoms = self._anoms_to_pandas(r_analysis)

        analysis = Analysis(EEG(self.tsubs(s, e)),
                            panoms=anoms[1], canoms=anoms[0])
        return (analysis)

    def analyze_stepwise(self, s, e, epoch=60):
        """Performs stepwise M-CAPA analysis upon a subset of the EEG
        with steps of size epoch.

        Args:
            s (float): Starting time
            e (float): Ending time
            epoch (int, optional): Size of each step in seconds. Defaults to 60.

        Returns:
            Analysis: An instance of class Analysis whose data is the
            subset of the EEG upon which M-CAPA analysis was performed,
            and whose canoms and panoms attributes are the collective
            and point anomalies found, respectively.
        """
        segm = self.tsubs(s, e)
        # Convert to R data frame
        with localconverter(ro.default_converter + pandas2ri.converter):
            r_eeg = ro.conversion.py2rpy(segm)

        # Create r eeg object
        r_eeg = ro.r.new("eeg", data=r_eeg)
        # Create r analysis object
        ro.r.stepwise(r_eeg, epoch, res=1, alpha=1)
        canoms = pd.read_csv("results/canoms.csv")
        panoms = pd.read_csv("results/panoms.csv")
        print("Done through reading")
        return (Analysis(EEG(segm), canoms, panoms))

    def draw(self, s, e, joint=False, width=1000, height=600):
        """Draws all EEG channels.

        Args:
             s (float): Starting time
            e (float): Ending time
            joint (bool, optional): Should EEG channels be plotted on top of one another? Defaults to False.
            width (int, optional): Width of the plot. Defaults to 1000.
            height (int, optional): Height of the plot. Defaults to 600.

        Returns:
            Figure (plotly): A plotly figure of all channels.
        """

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
