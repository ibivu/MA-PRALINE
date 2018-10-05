# Motif-Aware PRALINE (MA-PRALINE)

A motif-aware multiple sequence alignment tool. An article accompanying this software has been provisionally accepted for publication in PLOS Computational Biology; it can be found [here](https://github.com/ibivu/MA-PRALINE-extras/blob/master/article/article.pdf). 

# Installing MA-PRALINE

## Requirements

* Python 2.7 / Python 3.6 (earlier 3.x versions may also work, but have not been tested)
* An installation of PRALINE (https://github.com/ibivu/PRALINE)

## Instructions

You can install MA-PRALINE by cloning this repository and running (in a shell):

`python setup.py install`

A version will be uploaded to PyPI in the near future, which should make
MA-PRALINE installable with:

`pip install mapraline`

# Extras

There is an additional repository at https://github.com/ibivu/MA-PRALINE-extras containing
auxiliary code that was used to perform various analyses while developing MA-PRALINE.
Although no support is provided for this code, some people may find it useful to, e.g.
reparameterise MA-PRALINE or to calculate motif statistics.
