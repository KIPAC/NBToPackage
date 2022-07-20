# Software implementation, questions and answers


## 1. Having a separate class for every single plot seems like a bit of overkill, what are we gaining by doing things this way?

It gives us a way to configure the plots, and keep track of which configuration was used for each plot.  Basically it give us the provenance tracking & reproducibility.


## 2. The plot-specific function `_make_plot(self, cat)` function only take one argument, namely the input catalog, but the function the user sees `__call__(self, cat, **kwargs)` takes any number of argument through the `**kwargs` mechanism.  Why is this useful? 

This give us a way to latch and store the configuration used to make the plot.


## 3. What do we get from having the `PlotConfig` class?

It gives us a way to associated a name for a plot with a `Plotter` subclass and some configuration parameters.


## 4. What do we get from having the `PlotCollection` class?

It gives use a top-level thing to talk to, which keeps track of all the plots we have.


## 5. What did we come up with for dealing with common parameters?

We stuck them in a single file, to make it easy to find them.


## 6. What is your overall take on the final package?

Well, I'm pretty happy with it.

<!--  LocalWords:  kwargs kwargs
 -->
