# Original Notebooks, questions and answers

## 1. What is the notebook doing?

It is making a bunch of plots for a paper.  One really useful thing to notice is that all the plots are being made from the same input catalog.

## 2. How much redundancy / code duplication in the notebook that could be reduced?

There is a bit, but actually most of the plots are somewhat independent of each other.

## 3. What things are we are doing that make it less that ideal to be using a notebook?  Why?

A lot.

1. If we are going to be making plots for a paper, we really want things to be controlled and reproducible.  
2. This means we want to be doing things like keeping track of which data were used to make the plots, which version of the software was used to make the plots. 
3. Also, we are going to want to make it easy to standardize things like axis labels, axis limits, etc... 
4. In practice often just want to remake one plot, so we just re-run a couple cells out of order, it's cool that you can do this in a notebook, but it is really dangerous b/c it make things depend on the exact order you ran the cells.


## 4. How easy would it be to come back to this notebook in a year and expand on it?  Why?

It wouldn't be too bad, but the notebook is already pretty big, finding things and updating things is a bit of a pain.


## 5. How easy would it be to hand this notebook off to someone else?  Why?

I've seen worse.  The notebook is big and unwieldy, and there are a number of bits where it is implicitly picking up variables, but it has a pretty linear structure, so people can figure out what it is doing.
