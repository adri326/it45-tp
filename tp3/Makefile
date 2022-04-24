CC_FLAGS += -Wall -Wextra
LINK_FLAGS += -lm
ifdef DEBUG
CC_FLAGS += -g
else
CC_FLAGS += -O3
endif

CC_FLAGS += -g

ifdef VERBOSE
CC_FLAGS += -DVERBOSE
endif

BUILD_DIR = build
SRC_DIR = .
SRC = $(SRC_DIR)/little.c $(SRC_DIR)/main.c $(SRC_DIR)/test.c $(SRC_DIR)/time.c
OBJ = $(patsubst $(SRC_DIR)/%.c,$(BUILD_DIR)/%.o,$(SRC))

all: $(BUILD_DIR) $(BUILD_DIR)/little $(BUILD_DIR)/little-test $(BUILD_DIR)/little-time

clean:
	rm -f $(OBJ) $(BUILD_DIR)/test.o $(BUILD_DIR)/little $(BUILD_DIR)/little-test

format:
	clang-format -i little.c little.h test.c

report: $(BUILD_DIR)/report.pdf

.PHONY: all clean format report

$(BUILD_DIR)/report.pdf: report.tex
	pdflatex -interaction=nonstopmode -halt-on-error -shell-escape -output-directory=$(BUILD_DIR) report.tex

$(BUILD_DIR):
	mkdir -p $(BUILD_DIR)

$(BUILD_DIR)/little: $(BUILD_DIR)/main.o $(BUILD_DIR)/little.o
	$(CC) $(CC_FLAGS) $(LINK_FLAGS) $^ -o $@

$(BUILD_DIR)/little-test: $(BUILD_DIR)/test.o $(BUILD_DIR)/little.o
	$(CC) $(CC_FLAGS) $(LINK_FLAGS) -lcheck $^ -o $@

$(BUILD_DIR)/little-time: $(BUILD_DIR)/time.o $(BUILD_DIR)/little.o
	$(CC) $(CC_FLAGS) $(LINK_FLAGS) -lcheck $^ -o $@

$(OBJ): $(BUILD_DIR)/%.o: $(SRC_DIR)/%.c $(SRC_DIR)/little.h
	$(CC) $(CC_FLAGS) -c $< -o $@
