APP      = segment

CC ?= cc

SRCEXT   = c
SRCDIR   = src
OBJDIR   = obj
BINDIR   = bin

SRCS    := $(shell find . -name '*.$(SRCEXT)')
SRCDIRS := $(shell find . -name '*.$(SRCEXT)' -exec dirname {} \; | uniq)
OBJS    := $(patsubst %.$(SRCEXT),$(OBJDIR)/%.o,$(SRCS))

DEBUG    = -g
INCLUDES = -I./inc

# This 1990s source is K&R-flavoured; clang 15+ made several of its idioms
# hard errors.  None of them are bugs here (every implicitly-declared function
# in this tree returns int), so demote them back to warnings rather than
# rewrite the code.  No optimisation and no fast-math: the float arithmetic in
# region.c must not be reassociated or contracted.
WARNOFF  = -Wno-implicit-function-declaration \
           -Wno-implicit-int \
           -Wno-int-conversion \
           -Wno-deprecated-non-prototype \
           -Wno-extra-tokens \
           -Wno-macro-redefined \
           -Wno-return-type \
           -Wno-parentheses \
           -Wno-varargs \
           -Wno-format

# EXTRA_CFLAGS lets the flip() experiment be built without editing region.c:
#   make EXTRA_CFLAGS=-DFORCE_FLIP=0
CFLAGS   = -c -O0 -ffp-contract=off $(DEBUG) $(INCLUDES) $(WARNOFF) $(EXTRA_CFLAGS)
LDFLAGS  = -lm

.PHONY: all clean distclean

all: $(BINDIR)/$(APP)

$(BINDIR)/$(APP): buildrepo $(OBJS)
	@mkdir -p `dirname $@`
	@echo "Linking $@..."
	@$(CC) $(OBJS) $(LDFLAGS) -o $@

$(OBJDIR)/%.o: %.$(SRCEXT)
	@echo "Compiling $<..."
	@$(CC) $(CFLAGS) $< -o $@

clean:
	$(RM) -r $(OBJDIR)

distclean: clean
	$(RM) -r $(BINDIR)

buildrepo:
	@$(call make-repo)

define make-repo
   for dir in $(SRCDIRS); \
   do \
	mkdir -p $(OBJDIR)/$$dir; \
   done
endef
