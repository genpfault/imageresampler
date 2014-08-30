SOURCES_C   := $(wildcard *.c)
SOURCES_CPP := $(wildcard *.cpp)
OBJECTS_C   := $(patsubst %.c,bin/%.o,$(notdir $(SOURCES_C)))
OBJECTS_CPP := $(patsubst %.cpp,bin/%.o,$(notdir $(SOURCES_CPP)))
OBJECT_DIR = bin

bin/resampler: $(OBJECTS_C) $(OBJECTS_CPP)
	g++ $^ -o $@

$(OBJECT_DIR)/%.o : %.c
	$(CC) $(CPPFLAGS) $(CFLAGS) -c -o $@ $<

$(OBJECT_DIR)/%.o : %.cpp
	$(CC) $(CPPFLAGS) $(CFLAGS) -c -o $@ $<
