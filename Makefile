## This is Ningrui's thesis repo

current: target
-include target.mk
Ignore = target.mk

# -include makestuff/perl.def

vim_session:
	bash -cl "vmt"

######################################################################

## After_reviewing_papers.md Notes.md README.md todo.md
Sources += $(wildcard *.md)

######################################################################

## Code

Sources += $(wildcard Code/*.R Code/*.md)

## Depends on package
simulation.Rout: Code/simulation.R
	$(pipeR)

SeIR.Rout: Code/SeIR.R
	$(pipeR)

fitting.Rout: Code/fitting.R
	$(pipeR)

## foFits.md Notes about attempts to use fitode
foFits.Rout: Code/foFits.R
	$(pipeR)

######################################################################

tempfunc.Rout: Code/tempfunc.R
	$(pipeR)

hardtest_error.Rout: Code/hardtest_error.R tempfunc.rda
	$(pipeR)

######################################################################

## latex making example

subdirs += thesis figures

Sources += $(wildcard thesis/*.tex)

thesis/example.pdf: thesis/example.tex

alwayspipeR = defined

Sources += $(wildcard figures/*.R)
figures/test.Rout: figures/test.R

Ignore += subdirs

######################################################################

## JD stuff

multiFit.Rout: Code/multiFit.R
	$(pipeR)

oneFit.Rout: Code/oneFit.R
	$(pipeR)

compare.Rout: Code/compare.R
	$(pipeR)

## intervals.md
intervals.Rout: Code/intervals.R
	$(pipeR)

######################################################################

## JD older

simulate.Rout: jd/simulate.R
	$(pipeR)

functions.Rout: jd/functions.R
	$(pipeR)

######################################################################

## Modularize XNR's stuff? No progress 2023 Nov 09 (Thu)

peFuns.Rout: Code/peFuns.R

######################################################################

Ignore += .DS_Store

## Docs subdirectory
## docs/Makefile

pdirs += Proposal Midyear

## Outdated
Ignore += *.docs.pdf
%.docs.pdf: $(wildcard docs/*.*)
	cd docs && $(MAKE) $*.pdf
	$(CP) docs/$*.pdf $@

Ignore += $(pdirs)
alldirs += $(pdirs)

######################################################################

psync = $(pdirs:%=%.sync) 
psync: sync
	$(MAKE) $(psync)

######################################################################

### Makestuff

Sources += Makefile

Ignore += makestuff
msrepo = https://github.com/dushoff

Makefile: makestuff/00.stamp
makestuff/%.stamp:
	- $(RM) makestuff/*.stamp
	(cd makestuff && $(MAKE) pull) || git clone --depth 1 $(msrepo)/makestuff
	touch $@

-include makestuff/os.mk

-include makestuff/pipeR.mk
-include makestuff/texi.mk

-include makestuff/git.mk
-include makestuff/visual.mk
