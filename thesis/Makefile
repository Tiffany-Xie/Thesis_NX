## This is the thesis subdirectory of Thesis_NX

current: target
-include target.mk
Ignore = target.mk

# -include makestuff/perl.def

vim_session:
	bash -cl "vmt"

######################################################################

pardirs += figures
autopipeR = defined

Sources += Reference.bib
Sources += $(wildcard *.tex)

## draft.pdf: draft.tex

######################################################################

Sources += $(wildcard *.R)
Sources += $(wildcard figures/*.R)

tempfunc.Rout: ../Code/tempfunc.R
	$(pipeR)

figures/newtest.Rout: figures/test.R
	$(pipeR)

## It would be fun to write a pattern instead of repeating all of these rules
figures/difftypeFitting.Rout: figures/difftypeFitting.R tempfunc.rda #
	$(pipeR)

figures/sametypeFittingL.Rout: figures/sametypeFittingL.R tempfunc.rda #
	$(pipeR)

figures/sametypeFittingG.Rout: figures/sametypeFittingG.R tempfunc.rda #
	$(pipeR)

figures/parEst.Rout: figures/parEst.R tempfunc.rda #
	$(pipeR)

######################################################################

### Makestuff

Sources += Makefile

Ignore += makestuff
msrepo = https://github.com/dushoff

stufflink:
	ls ../makestuff/Makefile && ln -s ../makestuff
Makefile: makestuff/00.stamp
makestuff/%.stamp:
	- $(RM) makestuff/*.stamp
	(cd makestuff && $(MAKE) pull) || git clone $(msrepo)/makestuff
	touch $@

-include makestuff/os.mk

-include makestuff/texi.mk
-include makestuff/pdfpages.mk
-include makestuff/pipeR.mk

-include makestuff/git.mk
-include makestuff/visual.mk
