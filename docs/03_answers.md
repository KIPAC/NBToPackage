# Software design, questions and answers


## 1. What is the basic unit that we are going to work with?  A function that makes a single plot?  Several plots?

Probably a object that makes a single plot.  We want to keep in mind that we might want different version of the same plot, e.g., with
different limits, or a log version and a linear version.


## 2. How are we going to want to keep that of the plots we've made?  What bit of software knows about all the plots we've made?

Probably have a class with some internal dictionaries to keep track of things.  Probably want the dictionaries keyed by a name that we give to the plot, 
which could be the same as the filename it get written to.


## 3. How do we handle the balance between configurability and reproducibility?  How do we keep track of how a particular plot was made?

We should store the configuration when we make a plot.  Make have a base class that does that, and a derived class that makes the plot
from the configuration.


## 4. How do we manage common parameters?

Maybe just have a file for things like common labels, names, etc...


## 5. Do we want a command-line tool?  What options and arguments should it take?

Yes, it probably wants: path to input data, path to output locations, names of plots we want to make and a few other optional parameters.
