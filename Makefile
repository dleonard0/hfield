
CFLAGS += -g
CFLAGS += -O
CFLAGS += -Wall

ifdef NDEBUG
  CPPFLAGS += -DNDEBUG
endif

default: h

clean:
	rm -f h
