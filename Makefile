CC      ?= gcc
CFLAGS  ?= -O2 -std=gnu11 -Wall -Wextra -Wshadow -Wconversion -Wno-sign-conversion
CPPFLAGS = -Iinclude
LDLIBS   = -lm

BUILD := build
BIN   := $(BUILD)/fir_compare
SRC   := $(wildcard src/*.c)
OBJ   := $(patsubst src/%.c,$(BUILD)/%.o,$(SRC))
DEP   := $(OBJ:.o=.d)

.PHONY: all run test sweep clean

all: $(BIN)

$(BIN): $(OBJ)
	$(CC) $(CFLAGS) -o $@ $^ $(LDLIBS)

$(BUILD)/%.o: src/%.c | $(BUILD)
	$(CC) $(CFLAGS) $(CPPFLAGS) -MMD -MP -c -o $@ $<

$(BUILD):
	mkdir -p $(BUILD)

run: $(BIN)
	./$(BIN)

test: $(BIN)
	./$(BIN) --selftest

sweep: $(BIN)
	./$(BIN) --sweep

clean:
	rm -rf $(BUILD)

-include $(DEP)
