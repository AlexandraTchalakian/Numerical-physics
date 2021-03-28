CC = g++
CFLAGS = -Wall  -pedantic -ansi -g -std=c++11
EXEC_NAME = Gravitation
INCLUDES =
LIBS =
OBJ_FILES = PhysEngine.o Main.o
ZIP_FILES = $(addsuffix .cpp, $(basename $(OBJ_FILES))) $(addsuffix .hpp, $(basename $(OBJ_FILES))) configuration.in ConfigFile.hpp ConfigFile.tpp makefile GravitationScan
INSTALL_DIR = /usr/bin

all : $(EXEC_NAME)

clean :
	rm $(EXEC_NAME) $(OBJ_FILES)  

$(EXEC_NAME) : $(OBJ_FILES)
	$(CC) -o $(EXEC_NAME) $(OBJ_FILES) $(LIBS)

%.o: %.cpp
	$(CC) $(CFLAGS) $(INCLUDES) -o $@ -c $<

%.o: %.cc
	$(CC) $(CFLAGS) $(INCLUDES) -o $@ -c $<

%.o: %.c
	gcc $(CFLAGS) $(INCLUDES) -o $@ -c $<

install :
	cp $(EXEC_NAME) $(INSTALL_DIR)

zip :
	rm $(EXEC_NAME).zip
	mv PhysEngine.cpp PhysEngineSolution.cpp
	mv PhysEngineStudent.cpp PhysEngine.cpp
	zip $(EXEC_NAME).zip $(ZIP_FILES)
	mv PhysEngine.cpp PhysEngineStudent.cpp
	mv PhysEngineSolution.cpp PhysEngine.cpp

solution :
	rm -f $(EXEC_NAME)_solution.zip
	zip $(EXEC_NAME)_solution.zip $(ZIP_FILES)

