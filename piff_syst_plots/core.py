"""Core classes for piff_syst_plots"""

import os
from typing import Any
#import yaml

from collections import OrderedDict
from dataclasses import dataclass
import fitsio

from matplotlib.backends.backend_pdf import PdfPages


class Plotter:
    """Base class of object that make a single Figure"""

    plotters = OrderedDict()

    default_config = {}

    def __init_subclass__(cls, **kwargs):
        """
        Python 3.6+ provides a facility to automatically
        call a method (this one) whenever a new subclass
        is defined.  In this case we use that feature to keep
        track of all available pipeline stages, each of which is
        defined by a class.

        """
        super().__init_subclass__(**kwargs)
        cls.plotters[cls.__name__] = cls

    def __init__(self, **kwargs):
        self._config = self.default_config.copy()
        self._config.update(**kwargs)

    @property
    def config(self):
        """Configuration used to make plots"""
        return self._config

    def __call__(self, cat, **kwargs):
        """Main function, makes figure"""
        self._config.update(**kwargs)
        return self._make_plot(cat)

    def _make_plot(self, cat):
        """This function need to be implemented by the sub-class"""
        raise NotImplementedError()


@dataclass
class PlotConfig:
    """Simple class to carry around the configuration we use to make a plot"""

    plotname: str
    classname: str
    plot_config: dict[str, Any]

    def to_dict(self):
        """Return the PlotConfig as a dict"""
        return dict(
            plotname=self.plotname,
            classname=self.classname,
            plot_config=self.plot_config,
        )

    @classmethod
    def from_dict(cls, the_dict):
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


def add_pdf_metatada(pdffig, metadata):
    """Add provenance metadata to a pdf file"""
    pdf_metadata = pdffig.infodict()
    for key, value in metadata.items():
        pdf_metadata[key] = value


class PlotCollection:
    """A collection of plots we are making"""

    @staticmethod
    def get_version():
        """Return the version of the code used to make the plot"""
        return "v1"

    def __init__(self, data_url, output_dir, config=None):
        self._data_url = data_url
        self._catalog = fitsio.FITS(self._data_url)[1]
        self._catalog.upper = True
        self._catalog = self._catalog[:]
        self._output_dir = output_dir
        self._fig_dict = OrderedDict()
        self._plot_info_dict = OrderedDict()
        self._plotter_dict = OrderedDict()
        self.read_config(config)

    def read_config(self, config_url):
        """Get a list of PlotInfo objects from a yaml file"""
        if config_url is None:
            return
        with open(config_url, 'rt', encoding='utf-8') as config_file:
            #config_data = yaml.load(config_file)
            config_data = []
        for config in config_data:
            plot_info = PlotConfig.from_dict(config)
            self._plot_info_dict[plot_info.plotname] = plot_info

    def add_default_plotters(self):
        """Add a bunch of plotters"""
        for key in Plotter.plotters.keys():
            self._plot_info_dict[key] = PlotConfig(key, key, {})
            
    def make_plots(self, plot_names):
        """Make all the plots in the list"""
        if not plot_names:
            plot_names = self._plot_info_dict.keys()
        for plot_name in plot_names:
            plot_info = self._plot_info_dict[plot_name]
            the_plotter = plot_info.build_plotter()
            self._plotter_dict[plot_info.plotname] = the_plotter
            self._fig_dict[plot_info.plotname] = the_plotter(self._catalog)

    def _save_annotated(self, fig, fig_name, save_name):
        """Save the fig in a PDF with provencance metadata"""
        plot_info = self._plot_info_dict[fig_name]
        metadata = dict(
            version=self.get_version(),
            data_url=self._data_url,
            plot_function=plot_info.classname,
            plot_config=str(plot_info.plot_config),
        )
        fig.savefig(save_name)
        with PdfPages(save_name) as pdf:
            add_pdf_metatada(pdf, metadata)

    def show_fig(self, fig_name):
        """Show a single figure to the screen"""
        fig = self._fig_dict[fig_name]
        fig.show()

    def show_figs(self, fig_name_list=None):
        """Show all the figures to the screen"""
        if not fig_name_list:
            fig_name_list = self._fig_dict.keys()
        for fig_name in fig_name_list:
            self.show_fig(fig_name)

    def save_figs(self, fig_name_list=None, file_type="pdf", annotated=False):
        """Save all the figure in the list"""
        if not fig_name_list:
            fig_name_list = self._fig_dict.keys()

        try:
            os.makedirs(self._output_dir)
        except OSError:
            pass
        for fig_name in fig_name_list:
            fig = self._fig_dict[fig_name]
            if fig is None:
                return
            save_name = os.path.join(self._output_dir, f"{fig_name}.{file_type}")
            if file_type == "pdf" and annotated:
                self._save_annotated(fig, fig_name, save_name)
            else:
                fig.savefig(save_name)

    def add_plot(self, plot_name, plotter_class, **kwargs):
        """Add a plot to this collection"""
        plot_info = PlotConfig(plot_name, plotter_class.__name__, kwargs)
        the_plotter = plot_info.build_plotter()
        self._plot_info_dict[plot_info.plotname] = plot_info
        self._plotter_dict[plot_info.plotname] = the_plotter
        the_fig = the_plotter(self._catalog)
        self._fig_dict[plot_info.plotname] = the_fig
        return the_fig
