# Implementation of our package

I'm going to give a brief tour of the implementation I came up with.  Here are some question for you to consider:

1. Having a separate class for every single plot seems like a bit of overkill, what are we gaining by doing things this way?
2. The plot-specific function `_make_plot(self, cat)` function only take one argument, namely the input catalog, but the function the user sees `__call__(self, cat, **kwargs)` takes any number of argument through the `**kwargs` mechanism.  Why is this useful? 
3. What do we get from having the `PlotConfig` class?
4. What do we get from having the `PlotCollection` class?
5. What did we come up with for dealing with common parameters?
6. What is your overall take on the final package?

For what it's worth, here is [my take on those questions](04_answers.md)


<!--  LocalWords:  kwargs kwargs 04_answers.md
 -->
