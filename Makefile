
CFLAGS += -g
CFLAGS += -Wall

ifdef OPENMP
  CFLAGS += -O
  CFLAGS += -fopenmp
else
  CFLAGS += -Wno-unknown-pragmas
endif

ifdef NDEBUG
  CPPFLAGS += -DNDEBUG
endif

default: check

check: h t-h
	$(abspath t-h)
clean:
	rm -f h
