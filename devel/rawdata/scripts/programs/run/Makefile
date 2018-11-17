include ../*.mk

## check and recompile all *.o files upon which MAIN depends
$(ODIR)/%.o: $(SDIR)/%.cpp $(DEPS)
	$(CC) -c -o $@ $< $(AFLAGS) $(CFLAGS) &> $@.log

all: $(MAIN)
.PHONY: all

## link *.o files to create MAIN
$(MAIN): $(OBJ)
	$(CC) -o $@ $^ $(AFLAGS) $(LIBS) &> $@.log

## prevent creation of CLEAN file
.PHONY: cleanall cleanobj clean~ cleanprog

## clean directories of *.o and ~ files
cleanall: cleanobj clean~ cleanprog

cleanobj:
	-rm $(ODIR)/*.o

clean~:
	-rm *~

cleanprog:
	-rm $(MAIN) core