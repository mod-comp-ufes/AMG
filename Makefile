CC = gcc
CFLAGS = -O3
LDFLAGS = -lm
EXTENSION=.c
SRCDIR=src
OBJDIR=obj
SOURCES = $(shell cd $(SRCDIR) && find . -name "*$(EXTENSION)" | sed "s|^\./||")
OBJECTS = $(SOURCES:$(EXTENSION)=.o)
VPATH=$(dir $(addprefix $(SRCDIR)/,$(SOURCES)))
EXECUTABLE = program

all: $(EXECUTABLE)

debug:CFLAGS+=-g
debug: all

$(EXECUTABLE): $(addprefix $(OBJDIR)/,$(OBJECTS))
	$(CC) -o $@ $^ $(CFLAGS) $(LDFLAGS)

$(OBJDIR)/%.o: $(SRCDIR)/%$(EXTENSION)
	${shell mkdir -p ${dir $@}}
	$(CC) -c $(CFLAGS) $^ -o $@

.PHONY:clean
clean:
	rm -f $(addprefix $(OBJDIR)/,$(OBJECTS))