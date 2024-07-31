
2024 Jul 31 (Wed)
=================

[Here](https://github.com/eliminaterabies/R0paper/blob/main/public_data/intervals.rda) are some data you can try fitting as point data first; we will be curious how well they fit with dperlang, dlnorm, and dgamma. As we try things with real data, we should be making notes in a [stats journal](statsJournal.md) of everything we try, in order to maximize transparency and minimize bias.

Here are some pretty-complicated papers about thinking about intervals and data. You don't need to worry about them too much.
* [Charniga pre-print](https://arxiv.org/abs/2405.08841)
* [Park pre-print](https://www.medrxiv.org/content/10.1101/2024.01.12.24301247v1)

Ningrui is working on bounds for qperlang uniroot. A new idea would be to use qperlang with npe boxes and the slowest rate (for upper bound) or the highest rate: this is very non-conservative, but should work. The concern will be if it gives us ∞ more often than we like, but can try first.

In addition to trying that, we can also try to improve the current boundSearch method, which is working well for the parameters that Ningrui is using, but needs to be more flexible (and needs to not have magic numbers 🙂). Jonathan suggests a multiplicative approach with something like mean(exp(something)) where something is a vector symmetric around zero. To avoid ties, lower/upper bound can be multiplied by exp(∓ε).

Ningrui has been looking at the mpox data, and they are confusing. We could use point-based or interval-based estimation. Point-based estimation would involve just getting average dates; interval-based estimation is a little complicated if we have intervals at each end. It might be fun for Ningrui to do some math and write a function to describe the triangular distributions that result. We should also look at work by Abbott and Park and others where they may have code that we can leverage for that. Jonathan is thinking of finding some rabies data that may have a wider distribution and may be less problematic for point-based fitting. Ningrui will also look more at the COVID data.

Triangular distribution is something like Int_min^max of something about possible observation times. The probability of observing the data is then proportional to Int_min^max (pobs(int) * pperlang(int)).

2024 Jul 24 (Wed)
=================

inversion-based rperlang is now working.

We no longer need or or want rejection-based rperlang, but note that choosing c is hard, and probably requires some arbitrary limit (choose values where perlang is not too close to 1).

We talked about “normalizing” and “standardizing”. Standardizing in our field is the practice of giving things unit mean and variance (z = (x-μ)/σ). Normalizing is something that maybe we should be doing more: it means calculating a size statistic for a vector and dividing by it, to get a new vector whose size is 1. For densities, we would do something like d = d/(width*sum(d)), so that the new d would satisfy sum(width*d) = 1. It is always a good idea to: normalize vectors when appropriate; check whether the normalizing constant (the calculated size) is within expected limits.

Jonathan suggested some μ/κ-based upper bounds for the qperlang upper bound, these did not work!

Jonathan is supposed to look for data that we can fit.
* My Taiwan collaborator Andrei curated [some mpox data](https://github.com/aakhmetz/Mpox-IncubationPeriodSerialInterval-Meta2023/blob/main/SupplementaryFile1.xlsx). I will talk to him about it this week.
* Here is a [known set of covid interval data](https://zenodo.org/records/3940300)

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
