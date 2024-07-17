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

realfitting.Rout: Code/realfitting.R tempfunc.rda
	$(pipeR)

tempErrorPlot.Rout: Code/tempErrorPlot.R tempfunc.rda
	$(pipeR)


## Compare XNR's derived PDF to simulated density
testPEpdf.Rout: Code/testPEpdf.R
	$(pipeR)

## Working on fitting intervals 2024 Jul 03 (Wed)
Interval.Rout: Code/Interval.R tempfunc.rda
	$(pipeR)

######################################################################

Sources += $(wildcard data/*.csv)
data/WHO_covid.csv:
	wget -O $@ "https://covid19.who.int/WHO-COVID-19-global-data.csv"
WHO_covid.newdata:
	$(RM) data/WHO_covid.csv

WHO_covid.Rout: Code/WHO_covid.R data/WHO_covid.csv
	$(pipeR)

######################################################################

autopipeR = defined

tempfunc.Rout: Code/tempfunc.R
	$(pipeR)

hardtest_error.Rout: Code/hardtest_error.R tempfunc.rda
	$(pipeR)

## figures/mytest1.Rout: figures/mytest1.R

figures/mytest2.Rout: figures/mytest2.R tempfunc.rda
	$(pipeR)

figures/sametypeFittingL.Rout: figures/sametypeFittingL.R tempfunc.rda #
	$(pipeR)

figures/sametypeFittingR.Rout: figures/sametypeFittingR.R tempfunc.rda #
	$(pipeR)

figures/difftypeFitting.Rout: figures/difftypeFitting.R tempfunc.rda #
	$(pipeR)

figures/parEst.Rout: figures/parEst.R tempfunc.rda #
	$(pipeR)

######################################################################

## latex making example
Sources += $(wildcard *.tex)

subdirs += thesis figures

## Thesis is now made in thesis/ directory

autopipeR = defined

######################################################################

## Probably not making stuff here anymore 2024 Apr 05 (Fri)

Sources += $(wildcard figures/*.R)
figures/test.Rout: figures/test.R

Ignore += $(subdirs)

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

Ignore += docs
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
