
2024 Jul 17 (Wed)
=================

Ningrui has checked that naive (point-based) interval densities add up very close to 1; consider normalizing inside the code anyway.
* Or maybe not worry about it for now, since we will be moving on to real interval densities

You should be doing both of these things in many cases: check that it's close _and_ normalize anyway.

rperlang seems clever, but maybe not efficient. You could pick a bunch of parameters and calculate full densities to make a better choice for c.
c should be greater than max(dperlang(...)/dgamma(...) but not too large.

rperlang is rejection sampling, but there's another method that is often better, which involves inverting an estimated cdf.
* maybe try to build this with uniroot and pperlang

We should look at likelihoods and estimates of μ and κ, especially for gamma −> Erlang vs. gamma −> pseudoErlang

the density comparisons are cool; we might understand them a bit better if we look at some densities plotted on log scale

formula CDF is amazing, thanks Ningrui for not listening to Jonathan 🙂. We should talk about methods for inversion-based rperlang (although current rperlang seems good enough for now).

formula CDF will give us good interval densities!

Need to set or get histogram breaks and convert densities to counts using the right units – OR, there may be an option to make density-based histograms.
