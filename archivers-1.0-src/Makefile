#*************************************************************************
#  Sequential Online Archiving of Objective Vectors
#
VERSION=1.0
#                          Copyright (c) 2011
#         Manuel Lopez-Ibanez <manuel.lopez-ibanez@ulb.ac.be>
#             Joshua Knowles <j.knowles@manchester.ac.uk>
#                 Marco Laumanns <mlm@zurich.ibm.com>
#
# This program is free software (software libre); you can redistribute
# it and/or modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2 of the
# License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, you can obtain a copy of the GNU
# General Public License at http://www.gnu.org/copyleft/gpl.html
#
# IMPORTANT NOTE: Please be aware that the fact that this program is
# released as Free Software does not excuse you from scientific
# propriety, which obligates you to give appropriate credit! If you
# write a scientific paper describing research that made substantive use
# of this program, it is your obligation as a scientist to (a) mention
# the fashion in which this software was used in the Methods section;
# (b) mention the algorithm in the References section. The appropriate
# citation is:
# 
#     M. Lopez-Ibanez, J. Knowles, and M. Laumanns. On Sequential
#     Online Archiving of Objective Vectors. In R. Takahashi et al., 
#     editors, Evolutionary Multi-criterion Optimization (EMO 2011),
#     volume 6576 of Lecture Notes in Computer Science, pages 46-60.
#     Springer, Heidelberg, Germany, 2011.
# 
# Moreover, as a personal note, I would appreciate it if you would email
# manuel.lopez-ibanez@ulb.ac.be with citations of papers referencing this
# work so I can mention them to my funding agent and tenure committee.
#
# ----------------------------------------------------------------------
#
# Relevant literature:
#
# [1] M. Lopez-Ibanez, J. Knowles, and M. Laumanns. On Sequential
#     Online Archiving of Objective Vectors. In R. Takahashi et al., 
#     editors, Evolutionary Multi-criterion Optimization (EMO 2011),
#     volume 6576 of Lecture Notes in Computer Science, pages 46-60.
#     Springer, Heidelberg, Germany, 2011.
#
#*************************************************************************

DEBUG ?= 1

## Do we have svnversion?
ifeq ($(shell sh -c 'which svnversion 1> /dev/null 2>&1 && echo y'),y)
  ## Is this a working copy?
  ifneq ($(shell sh -c 'LC_ALL=C svnversion -n .'),exported)
    $(shell sh -c 'svnversion -n . > svn_version')
  endif
endif

## Set version information:
SVN_REV = $(shell sh -c 'cat svn_version 2> /dev/null')

MAIN = archiver
SRCS = archiver.cpp Archive.cpp Solution.cpp rng.c mt.c
HDRS = *.h
OBJS = $(patsubst %.cpp,%.o,$(patsubst %.c,%.o,$(SRCS)))

DIST_SRC_FILES = Makefile Makefile.lib \
		 README LICENSE PISA_LICENSE.txt \
		 svn_version \
		 $(SRCS) $(HDRS) \
		 $(HV_SRCS) $(HV_HDRS)

DIST_SRC       = archivers-$(VERSION)-src

################################################################################
## Configure the compiler / linker:

## Global list of CPP flags
CPPFLAGS = -D DEBUG=$(DEBUG) -D VERSION='"$(VERSION)"'

ifneq ($(DEBUG), 0)
DEBUGFLAGS = -g3
CPPFLAGS += -DMALLOC_CHECK_=3
endif

CXXFLAGS = $(DEBUGFLAGS) -Wall -Wextra
CFLAGS = $(CXXFLAGS)

.PHONY: all clean dist test default

all: clean $(MAIN)

clean:
	$(call ECHO,---> Removing $(MAIN) <---)
	@$(RM) $(MAIN)
	$(call ECHO,---> Removing object files <---)
	@$(RM) $(OBJS) $(HV_OBJS)
	$(call ECHO,---> Removing $(HV_LIB) <---)
	@$(RM) $(HV_LIB)
	$(call ECHO,---> Removing backup files <---)
	@$(RM) *~

dist : all
	@(rm -f ../$(DIST_SRC).tar.gz && mkdir -p ../$(DIST_SRC) \
	&& rsync -rlpC --relative --exclude=.svn $(DIST_SRC_FILES) ../$(DIST_SRC)/ \
	&& cd .. \
	&& tar cf - $(DIST_SRC) | gzip -f9 > $(DIST_SRC).tar.gz \
	&& rm -rf ./$(DIST_SRC)/* && rmdir ./$(DIST_SRC)/ \
	&& echo "$(DIST_SRC).tar.gz created." && cd $(PWD) )

#----------------------------------------------------------------------
# Include actual HV code:
include Makefile.lib

$(MAIN): $(OBJS) $(HV_LIB)
	$(call ECHO,---> Building $@ version $(VERSION) <---)
	${CXX} ${CXXFLAGS} $^ -o $@

#----------------------------------------------------------------------
# Dependencies:
archiver.o: archiver.cpp eApproxArchive.h eParetoArchive.h SPEA2Archive.h NSGA2Archive.h HVArchive.h AdaptiveGridArchive.h MultilevelGridArchive.h Random.h 

Archive.o: Archive.cpp Archive.h

Solution.o: Solution.cpp Solution.h

mt.o: mt.c gsl_rng.h gsl_types.h

rng.o : rng.c gsl_rng.h gsl_types.h


