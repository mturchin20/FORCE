MKFILE_PATH := $(abspath $(lastword $(MAKEFILE_LIST)))
CURRENT_PATH := $(dir $(abspath $(MAKEFILE_LIST)))
CURRENT_DIR := $(notdir $(patsubst %/,%,$(dir $(MKFILE_PATH))))
HTML_FILES := $(patsubst %.Rmd, %.html ,$(wildcard $(CURRENT_PATH)analysis/*.Rmd)) \
              $(patsubst %.md, %.html ,$(wildcard $(CURRENT_PATH)analysis/*.md))
#HTML_FILES_OUTPUT := $(patsubst /analysis/, /docs/ , $(HTML_FILES))

all: html

html: $(HTML_FILES)

%.html: %.Rmd
	R --slave -e "set.seed(100);rmarkdown::render('$<')"
	cp -p $(CURRENT_PATH)analysis/$(notdir $@) $(CURRENT_PATH)docs/$(notdir $@) 

%.html: %.md
	R --slave -e "set.seed(100);rmarkdown::render('$<')"
	cp -p $(CURRENT_PATH)analysis/$(notdir $@) $(CURRENT_PATH)docs/$(notdir $@) 

.PHONY: clean
clean:
	$(RM) $(patsubst /analysis/, /docs/ , $(HTML_FILES)) 

