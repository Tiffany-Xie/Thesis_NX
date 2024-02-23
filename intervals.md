Park…Abbott Preprint (for reference)
* https://www.medrxiv.org/content/10.1101/2024.01.12.24301247v1

Our goal will be to fit to interval data.

We will do this under two paradigms: the “simple rounding paradigm”, where we assume that 2d means between 1.5d and 2.5d, and the “start/end paradigm”, where we take into account possible starting and ending times. We may never get to fitting the second paradigm, but we should certainly simulate from it.

We are going to use different distributions to generate data: probably lognormal and gamma.

We will fit things in three ways: gamma, Erlang, and pseudo-Erlang. To do the Erlang we can just fit the two shape parameters that flank the gamma-estimated shape parameter.

We could also compare naive (dgamma-based) gamma fitting to fancy (pgamma-based gamma fitting).
