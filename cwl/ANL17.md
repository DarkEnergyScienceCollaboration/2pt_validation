## Sprint at ANL 2017 

What we want to do is to convert the [Cookbook](https://github.com/LSSTDESC/2pt_validation/blob/master/doc/Cookbook.md)
instructions into a runnable CWL pipeline.

There are a number of issues:

* The process has evolved into something pretty messy over the past
  year. Hopefully this pipeline might put things back in order

* We didn't propagate documentatin properly. For example, there was a makefile for `CoLoRe` in 
NERSC directory that hasn't been updated over a year. We need to prevent every new person rediscovering how to build CoLoRe. 

* Do we want to have "officially" compiled versions of all codes? Do we want to containerize them? What is the right place to put them.

* How does one run the pipeline run part-way? Say I already have CoLore Outputs and just want to run from fastcat onwards?

* What is the best way to specify configuration parameters? Are they flat-files or CWL objects?

* How do we run a set of experiments where we vary some parameters and have a clear idea what run and what files were produced where?

