"""Core classes for piff_syst_plots"""

from __future__ import annotations

import os
from typing import Any, Optional

# import yaml

from collections import OrderedDict
from dataclasses import dataclass
import fitsio

import matplotlib as mpl
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.figure import Figure


def set_interactive():
    mpl.use('qtagg')


def set_paper_style():
    # set some plotting defaults
    mpl.rc(('lines', 'axes') , linewidth=2)
    mpl.rc(('xtick', 'ytick'), labelsize=15)
    mpl.rc(('xtick.major', 'ytick.major'), width=2)
    mpl.rcParams['axes.labelsize'] = 18
    mpl.rcParams['font.family'] = 'serif'
    mpl.rcParams['mathtext.fontset'] = 'stix'
    mpl.rcParams['legend.fontsize'] = 15
    mpl.rcParams['font.size'] = 18


def set_talk_style():
    # set some plotting defaults
    mpl.rc(('lines', 'axes') , linewidth=2)
    mpl.rc(('xtick', 'ytick'), labelsize=15)
    mpl.rc(('xtick.major', 'ytick.major'), width=2)
    mpl.rcParams['axes.labelsize'] = 26
    mpl.rcParams['font.family'] = 'serif'
    mpl.rcParams['mathtext.fontset'] = 'stix'
    mpl.rcParams['legend.fontsize'] = 18
    mpl.rcParams['font.size'] = 26


STYLES = dict(
    paper=set_paper_style,
    talk=set_talk_style,
)


def set_style(style):
    the_func = STYLES[style]
    the_func()


class Plotter:
    """Base class of object that make a single Figure

    This class really just keeps track of the options that
    you use to make the plot, via the `self._config` dictionary.

    Derived classes should override the `_make_plot(self, cat)`
    function to make and return a `matplotlib.figure.Figure`
    """

    plotters = OrderedDict()

    default_config = {}

    def __init_subclass__(cls, **kwargs):
        """Python 3.6+ provides a facility to automatically
        call a method (this one) whenever a new subclass
        is defined.  In this case we use that feature to keep
        track of all available plots, each of which is
        defined by a class.

        """
        super().__init_subclass__(**kwargs)
        cls.plotters[cls.__name__] = cls

    def __init__(self, **kwargs):
        """Constructor, just copy kwargs to self._config"""
        self._config = self.default_config.copy()
        self._config.update(**kwargs)

    @property
    def config(self) -> dict[str, Any]:
        """Configuration used to make plots"""
        return self._config

    def __call__(self, cat, **kwargs) -> Figure:
        """Main function, makes figure"""
        self._config.update(**kwargs)
        return self._make_plot(cat)

    def _make_plot(self, cat) -> Figure:
        """This function need to be implemented by the sub-class"""
        raise NotImplementedError()


@dataclass
class PlotConfig:
    """Simple class to carry around the configuration we use to make a plot

    This includes:

    Parameters
    ----------
    plotname : str
        The name of the plot, used to refer to it, and also used to make
        the output file name

    classname : str
        The name of the python class that will actually make the plot

    plot_config : dict[str, Any]
        The configuration options for making the plot
    """

    plotname: str
    classname: str
    plot_config: dict[str, Any]

    def to_dict(self) -> dict[str, Any]:
        """Return the PlotConfig as a dict"""
        return dict(
            plotname=self.plotname,
            classname=self.classname,
            plot_config=self.plot_config,
        )

    @classmethod
    def from_dict(cls, the_dict) -> PlotConfig:
        """Create a PlotConfig object from a dict"""
        return cls(
            plotname=the_dict["plotname"],
            classname=the_dict["classname"],
            plot_config=the_dict["plot_config"],
        )

    def build_plotter(self) -> Plotter:
        """Return a configured Plotter"""
        plotter_class = Plotter.plotters.get(self.classname)
        if plotter_class is None:
            raise KeyError(f"Could not find a plotter class called {self.classname}")
        return plotter_class(**self.plot_config)


def add_pdf_metatada(pdffig: PdfPages, metadata: dict[str, str]) -> None:
    """Add provenance metadata to a pdf file"""
    pdf_metadata = pdffig.infodict()
    for key, value in metadata.items():
        pdf_metadata[key] = value


class PlotCollection:
    """A collection of plots we are making

    Individual plots are refered to by `plot_name`, which
    is used as the key to all the internal dictionaries.
    """

    @staticmethod
    def get_version() -> str:
        """Return the version of the code used to make the plot"""
        return "v1"

    def __init__(self, data_url: str, output_dir: str, config: Optional[str] = None):
        """Constructor

        Parameters
        ----------
        data_url : str
            The path to the data being used to make these plots

        output_dir: str
            The path to where we want to store these plots

        config : str
            The path to a yaml file we can use to (re)-configure the plots
        """
        self._data_url = data_url
        self._catalog = fitsio.FITS(self._data_url)[1]
        self._catalog.upper = True
        self._catalog = self._catalog[:]
        self._output_dir = output_dir
        self._fig_dict = OrderedDict()
        self._plot_info_dict = OrderedDict()
        self._plotter_dict = OrderedDict()
        self.read_config(config)

    def read_config(self, config_url: str) -> None:
        """Get a list of PlotInfo objects from a yaml file"""
        if config_url is None:
            return
        with open(config_url, "rt", encoding="utf-8") as config_file:
            # config_data = yaml.load(config_file)
            config_data = []
        for config in config_data:
            plot_info = PlotConfig.from_dict(config)
            self._plot_info_dict[plot_info.plotname] = plot_info

    def add_default_plotters(self) -> None:
        """Add a bunch of plotters"""
        for key in Plotter.plotters.keys():
            self._plot_info_dict[key] = PlotConfig(key, key, {})

    def make_plots(self, plot_names: str) -> None:
        """Make all the plots in the list"""
        if not plot_names:
            plot_names = self._plot_info_dict.keys()
        for plot_name in plot_names:
            plot_info = self._plot_info_dict[plot_name]
            the_plotter = plot_info.build_plotter()
            self._plotter_dict[plot_info.plotname] = the_plotter
            self._fig_dict[plot_info.plotname] = the_plotter(self._catalog)

    def _save_annotated(self, fig: Figure, plot_name: str, save_name: str) -> None:
        """Save the fig in a PDF with provencance metadata"""
        plot_info = self._plot_info_dict[plot_name]
        metadata = dict(
            version=self.get_version(),
            data_url=self._data_url,
            plot_function=plot_info.classname,
            plot_config=str(plot_info.plot_config),
        )
        fig.savefig(save_name)
        with PdfPages(save_name) as pdf:
            add_pdf_metatada(pdf, metadata)

    def show_fig(self, plot_name: str) -> None:
        """Show a single figure to the screen"""
        fig = self._fig_dict[plot_name]
        fig.show()

    def show_figs(self, plot_name_list: Optional[list[str]] = None):
        """Show all the figures to the screen"""
        if not plot_name_list:
            plot_name_list = self._fig_dict.keys()
        for plot_name in plot_name_list:
            self.show_fig(plot_name)
        mpl.pyplot.show()

    def save_figs(
        self,
        plot_name_list: Optional[list[str]] = None,
        file_type: str = "pdf",
        annotated: bool = False,
    ) -> None:
        """Save all the figures in the list"""
        if not plot_name_list:
            plot_name_list = self._fig_dict.keys()

        try:
            os.makedirs(self._output_dir)
        except OSError:
            pass
        for plot_name in plot_name_list:
            fig = self._fig_dict[plot_name]
            if fig is None:
                return
            save_name = os.path.join(self._output_dir, f"{plot_name}.{file_type}")
            if file_type == "pdf" and annotated:
                self._save_annotated(fig, plot_name, save_name)
            else:
                fig.savefig(save_name)

    def add_plot(self, plot_name: str, plotter_class: type, **kwargs) -> Figure:
        """Add a plot to this collection

        Parameters
        ----------
        plot_name : str
            Name we will use for the plot
            this is also used to make the output file name

        plotter_class : type
            Class actually used to make the plot

        Keywords are passed to the plotter_class to configure it
        """
        plot_info = PlotConfig(plot_name, plotter_class.__name__, kwargs)
        the_plotter = plot_info.build_plotter()
        self._plot_info_dict[plot_info.plotname] = plot_info
        self._plotter_dict[plot_info.plotname] = the_plotter
        the_fig = the_plotter(self._catalog)
        self._fig_dict[plot_info.plotname] = the_fig
        return the_fig
