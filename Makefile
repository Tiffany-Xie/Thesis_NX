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

Sources += $(wildcard Code/*.R)

## Depends on package
simulation.Rout: Code/simulation.R
	$(pipeR)

SeIR.Rout: Code/SeIR.R
	$(pipeR)

fitting.Rout: Code/fitting.R
	$(pipeR)

######################################################################

## JD stuff

multiFit.Rout: Code/multiFit.R
	$(pipeR)

oneFit.Rout: Code/oneFit.R
	$(pipeR)

compare.Rout: Code/compare.R
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

pdirs += docs

## Proposal.docs.pdf: docs/Proposal.tex
Ignore += *.docs.pdf
%.docs.pdf: $(wildcard docs/*.*)
	cd docs && $(MAKE) $*.pdf
	$(CP) docs/$*.pdf $@

Ignore += $(pdirs)
alldirs += $(pdirs)

######################################################################

propsub: docs/Proposal.tex.16571ad84.oldfile

## Proposal.ld.docs.pdf: docs/Proposal.ld.tex
## docs/Proposal.ld.tex: docs/Proposal.tex

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

-include makestuff/git.mk
-include makestuff/visual.mk
