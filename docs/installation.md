# installation

### Installing the package

Let's assume that you have installed conda or miniconda or mamba.  We are going to
use that to create an environment for this package.  We will add jupyter and fitsio
using conda because that seems to work better, and add this package using pip, because
that works better. 

	conda create -n nb_to_pack -c conda-forge -y python=3.10.5 jupyter fitsio
	conda activate nb_to_pack
	git clone https://github.com/KIPAC/NBToPackage.git
	cd NBToPackage
	pip install -e .


### Running notebooks

You can start up a jupyter notebook session (pointed at the directory with the notebooks) like so:

	jupyter-notebok nb &

From there can launch either notebook.  Note that we saved a version of the old notebook that was run, while we saved a 'clean' version of the new-notebook, to make it clearer what the new code looks like.

### Running the command line tool

Ok, let's try make the plots from the command line tool we build to replace the notebook.
First let's see how to run the tool.

	bin/piff-syst-plots -h

	usage: piff_syst_plots [-h] -d FILE [-o DIR] [-c FILE] [-f] [-a] [-l]
		[NAMES ...]

	positional arguments:
		NAMES                 Names of plots to draw, leave blank for all plots

	options:
		-h, --help            show this help message and exit
		-d FILE, --data FILE  Path to data
		-o DIR, --output DIR  Output directory
		-c FILE, --config FILE Yaml configuration file
		-p , --plot_type      Plot file type: ['pdf', 'png', 'eps', 'jpg']
		-a, --annotate        Store provenance metdata with pdf
		-l, --list_plotters   list plotters


Ok, let's make all the plots using the `data/piff_sample_cat.fits` file provided.

	bin/piff_syst_plots -d data/piff_sample_cat.fits -o plots -p pdf
	






	


