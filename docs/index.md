# Moving from jupyter notebook to a python package

We are going to be spending about and hour and a half working through how we might design code 
as a python package rather than as a jupyter notebook.

We are going to be working from a concrete example, but, given the time available, we will have to 
focus on big-picture points.

## Objectives and content

The objectives here are pretty simple
 
  1. understand the advantages moving from a jupyter-notebook analysis to a python package,
  2. understand how to design a really nice python package,
  3. understand how to build a really nice python package, a to get a sense of the tools available to help with that.

To reach these objectives we are going to break this down as follows

  1. start with a pretty typical notebook that makes a bunch of plots for a paper, and characterize what the notebook is doing,
  2. rethink what we want from the software, 
  3. come up with some design that gets us more of what we want,
  4. look at an example implementation of a python package that replaces the original notebook and see how it address what we came up with in 2,3,4 and 5.

## Pre-course materials

Please have a look at the slides from Dark Energy School [session on design code libraries and abstract public interfaces (APIs)](https://lsstdesc.org/pages/DESchool.html#JarvisZuntz), [slides are here](https://lsstdesc.org/assets/pdf/docs/ParisDESchoolAPIDesign.pdf)

If you want to set up the software package and try running both versions
of the notebook and the command line tool for comparison, here are
instructions on how to do so.   It should only take 10-15' minutes following these instructions:

[Installation instructions](installation.md)

On the other hand, you can also just follow along this documentation:

[https://kipac.github.io/NBToPackage/](https://kipac.github.io/NBToPackage/)

## Course Contents

First off, you probably want a quick round of introductions with the
people you are sitting or sharing a breakout room with.

Here is a handy link to a google doc template with the questions
copied out, you can do "File -> make a copy" to make your own doc that
you can use to work together to come up with answers to some of those
questions:
[doc template](https://docs.google.com/document/d/1mV_T4pvxMIC9lSOCJpREhpsWqbZGWYnC1OsEHrrL3hA/edit?usp=sharing)
I'd recommend making one copy of the file, and setting it up so you can
all edit it.

Ok, Let's get started, here are some times for guidance:

25' Total:  10' Intro/Overhead, 10' Breakout discussion, 5' Full group,
[Original Notebook](01_starting_with_the_original.md)

30' Total:  5' Intro/Overhead, 15' Breakout discussion, 10' Full group
[Software Goals](02_goals_for_our_software.md)
[Designing our Package](03_designing_our_package.md)


35' Total: 10' Intro/Overhead 15' Breakout discussion, 10' Full group
[Implementation of our Package](04_implementation_of_our_package.md)

Extras: [Python package tooling](extra_python_tooling.md)



<!--  LocalWords:  jupyter JarvisZuntz 02_goals_for_our_software.md
 -->
<!--  LocalWords:  01_starting_with_the_original.md
 -->
<!--  LocalWords:  03_designing_our_package.md
 -->
<!--  LocalWords:  04_implementation_of_our_package.md
 -->
<!--  LocalWords:  extra_python_tooling.md
 -->
