CC = gcc
#For older gcc, use -O3 or -O2 instead of -Ofast
# CFLAGS = -lm -pthread -Ofast -march=native -funroll-loops -Wno-unused-result

# Use -Ofast with caution. It speeds up training, but the checks for NaN will not work
# (-Ofast turns on --fast-math, which turns on -ffinite-math-only,
# which assumes everything is NOT NaN or +-Inf, so checks for NaN always return false
# see https://gcc.gnu.org/onlinedocs/gcc/Optimize-Options.html)
# CFLAGS = -lm -pthread -Ofast -march=native -funroll-loops -Wall -Wextra -Wpedantic

CFLAGS = -lm -pthread -O3 -march=native -funroll-loops -Wall -Wextra -Wpedantic
BUILDDIR := build
SRCDIR := src
OBJDIR := $(BUILDDIR)

OBJ := $(OBJDIR)/vocab_count.o $(OBJDIR)/cooccur.o $(OBJDIR)/shuffle.o $(OBJDIR)/glove.o $(OBJDIR)/glove_momentum_opt.o $(OBJDIR)/motiva.o $(OBJDIR)/cooccur_persuser.o $(OBJDIR)/mergeuserdata.o $(OBJDIR)/per_file_corr.o  $(OBJDIR)/realmerge.o $(OBJDIR)/shuffcipdata.o $(OBJDIR)/trainincip.o $(OBJDIR)/traincip_triple.o $(OBJDIR)/ciplogandf
HEADERS := $(SRCDIR)/common.h
MODULES := $(BUILDDIR)/vocab_count $(BUILDDIR)/cooccur $(BUILDDIR)/shuffle $(BUILDDIR)/glove $(BUILDDIR)/glove_momentum_opt $(BUILDDIR)/cooccur_persuser $(BUILDDIR)/mergeuserdata $(BUILDDIR)/per_file_corr $(BUILDDIR)/realmerge $(BUILDDIR)/shuffcipdata $(BUILDDIR)/trainincip $(BUILDDIR)/traincip_triple $(BUILDDIR)/ciplogandf


all: dir $(OBJ) $(MODULES)
dir :
	mkdir -p $(BUILDDIR)
$(BUILDDIR)/glove : $(OBJDIR)/glove.o $(OBJDIR)/common.o
	$(CC) $^ -o $@ $(CFLAGS)
$(BUILDDIR)/glove_momentum_opt : $(OBJDIR)/glove_momentum_opt.o $(OBJDIR)/common.o
	$(CC) $^ -o $@ $(CFLAGS)
$(BUILDDIR)/motiva : $(OBJDIR)/motiva.o $(OBJDIR)/common.o
	$(CC) $^ -o $@ $(CFLAGS)
$(BUILDDIR)/shuffle : $(OBJDIR)/shuffle.o $(OBJDIR)/common.o
	$(CC) $^ -o $@ $(CFLAGS)
$(BUILDDIR)/cooccur : $(OBJDIR)/cooccur.o $(OBJDIR)/common.o
	$(CC) $^ -o $@ $(CFLAGS)

$(BUILDDIR)/trainincip : $(OBJDIR)/trainincip.o $(OBJDIR)/common.o
	$(CC) $^ -o $@ $(CFLAGS)

$(BUILDDIR)/traincip_triple : $(OBJDIR)/traincip_triple.o $(OBJDIR)/common.o
	$(CC) $^ -o $@ $(CFLAGS)

$(BUILDDIR)/realmerge : $(OBJDIR)/realmerge.o $(OBJDIR)/common.o
	$(CC) $^ -o $@ $(CFLAGS)

$(BUILDDIR)/mergeuserdata : $(OBJDIR)/mergeuserdata.o $(OBJDIR)/common.o
	$(CC) $^ -o $@ $(CFLAGS)

$(BUILDDIR)/per_file_corr : $(OBJDIR)/per_file_corr.o $(OBJDIR)/common.o
	$(CC) $^ -o $@ $(CFLAGS)

$(BUILDDIR)/shuffcipdata : $(OBJDIR)/shuffcipdata.o $(OBJDIR)/common.o
	$(CC) $^ -o $@ $(CFLAGS)
 
$(BUILDDIR)/cooccur_persuser : $(OBJDIR)/cooccur_persuser.o $(OBJDIR)/common.o
	$(CC) $^ -o $@ $(CFLAGS)

$(BUILDDIR)/ciplogandf : $(OBJDIR)/ciplogandf.o $(OBJDIR)/common.o
	$(CC) $^ -o $@ $(CFLAGS)
	
$(BUILDDIR)/vocab_count : $(OBJDIR)/vocab_count.o $(OBJDIR)/common.o
	$(CC) $^ -o $@ $(CFLAGS)
$(OBJDIR)/%.o : $(SRCDIR)/%.c $(HEADERS)
	$(CC) -c $< -o $@ $(CFLAGS)
.PHONY: clean
clean:
	rm -rf $(BUILDDIR)
