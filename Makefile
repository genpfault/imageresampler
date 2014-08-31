CPPFLAGS = -Wall -Wextra
SOURCES_CPP := $(wildcard *.cpp)
OBJECTS_CPP := $(patsubst %.cpp,bin/%.o,$(notdir $(SOURCES_CPP)))
OBJECT_DIR = bin

bin/resampler: $(OBJECTS_CPP)
	g++ $^ -o $@

$(OBJECT_DIR)/%.o : %.cpp
	$(CC) $(CPPFLAGS) $(CFLAGS) -c -o $@ $<

clean:
	rm -f bin/resampler $(OBJECT_DIR)/*.o
