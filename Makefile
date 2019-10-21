
ifeq ($(strip $(CC)),)
CC := gcc
endif
CC += -Wall -MP -MD 
INCMLFLG := -D INC_ML_CMPL
TESTFLG := -D TEST
DIR := src/c

SPECCMPL := 

SOURCES := $(wildcard $(DIR)/*.c)
OBJECTS := $(patsubst $(DIR)/%.c,  $(DIR)/%.o, $(SOURCES))

.PHONY: dactal
# test_tree: SPECCMPL = -D TEST_TREE
# test_rf2: SPECCMPL = -D TEST_RF2 

dactal: $(OBJECTS)
	$(CC) $(SPECCMPL) $^ -lm -o $@

$(DIR)/%.o: $(DIR)/%.c
	$(CC) $(SPECCMPL) -I$(DIR) -c $< -o $@

clean:
	rm -f $(DIR)/*.o 
	rm -f $(DIR)/*.d
	rm -f dactal
