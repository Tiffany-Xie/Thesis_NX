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
## example.pdf: example.tex

## draft.pdf: draft.tex

Sources += $(wildcard *.R)
## figures/test.Rout: figures/test.R
## figures/excellent.Rout: figures/excellent.R
figures/newtest.Rout: figures/test.R
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
