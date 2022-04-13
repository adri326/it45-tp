CC_FLAGS += -Wall -Wextra
LINK_FLAGS += -lm
ifdef DEBUG
CC_FLAGS += -g
else
CC_FLAGS += -O3
endif

ifdef VERBOSE
CC_FLAGS += -DVERBOSE
endif

BUILD_DIR = build
SRC_DIR = .
SRC = $(SRC_DIR)/little.c
OBJ = $(patsubst $(SRC_DIR)/%.c,$(BUILD_DIR)/%.o,$(SRC))
OBJ_TEST = $(patsubst $(SRC_DIR)/%.c,$(BUILD_DIR)/%-test.o,$(SRC))

all: $(BUILD_DIR) $(BUILD_DIR)/little $(BUILD_DIR)/little-test

clean:
	rm -f $(OBJ) $(OBJ_TEST) $(BUILD_DIR)/test.o $(BUILD_DIR)/little $(BUILD_DIR)/little-test

format:
	clang-format -i little.c little.h test.c

$(BUILD_DIR):
	mkdir -p $(BUILD_DIR)

$(BUILD_DIR)/little: $(BUILD_DIR) $(OBJ)
	$(CC) $(CC_FLAGS) $(LINK_FLAGS) $(OBJ) -o $@

$(BUILD_DIR)/little-test: $(BUILD_DIR) $(OBJ_TEST) $(BUILD_DIR)/test.o
	$(CC) $(CC_FLAGS) $(LINK_FLAGS) -lcheck $(OBJ_TEST) $(BUILD_DIR)/test.o -o $@

$(OBJ): $(BUILD_DIR)/%.o: $(SRC_DIR)/%.c
	$(CC) $(CC_FLAGS) -c $< -o $@

$(BUILD_DIR)/test.o: $(SRC_DIR)/test.c
	$(CC) $(CC_FLAGS) -DUNIT_TEST -c $< -o $@

$(OBJ_TEST): $(BUILD_DIR)/%-test.o: $(SRC_DIR)/%.c
	$(CC) $(CC_FLAGS) -DUNIT_TEST -c $< -o $@
